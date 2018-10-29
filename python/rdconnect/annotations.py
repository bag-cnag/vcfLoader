from rdconnect import utils, expr

def importVCF(hc, sourcePath, destinationPath, nPartitions):
    """ Imports input vcf and annotates it with general annotations (samples, freqInt, pos, alt, ref)
          :param HailContext hc: The Hail context
          :param String sourcePath: Annotation table path
          :param String destinationPath: Path where the loaded annotation table will be put
          :param String nPartitions: Number of partitions
    """
    try:
        print ("reading vcf from "+ sourcePath)
        vcf = hc.import_vcf(str(sourcePath),force_bgz=True,min_partitions=nPartitions).split_multi()
        print ("writing vds to" + destinationPath)
        vcf.annotate_variants_expr(expr.annotationsVariants()) \
           .annotate_variants_expr(expr.annotationsFreqInt()) \
           .write(destinationPath,overwrite=True)
        return True
    except ValueError:
        print (ValueError)
        return "Error in importing vcf"

def importSomatic(hl, file_paths, destination_path, num_partitions):
    nFiles = len(file_paths)
    if(nFiles > 0) :
        try:
            merged = hl.split_multi(hl.import_vcf(file_paths[0],force_bgz=True,min_partitions=num_partitions))
            merged = annotateSomatic(hl,merged)
            for file_path in file_paths[1:]:
                print("File path -> " + file_path)
                dataset = hl.split_multi(hl.import_vcf(file_path,force_bgz=True,min_partitions=num_partitions))
                dataset = annotateSomatic(hl,dataset)
                merged = mergeSomatic(merged,dataset)
            merged.write(destination_path,overwrite=True)
        except ValueError:
            print("Error in loading vcf")
    else:
        print("Empty file list")

def mergeSomatic(dataset, other):
    tdataset = dataset.rows()
    tdataset.show()
    tother = other.rows()
    tother.show()
    joined = tdataset.join(tother,"outer")
    return joined.transmute(samples=joined.samples.union(joined.samples_1))
    
    
def importDbNSFPTable(hl, sourcePath, destinationPath, nPartitions):
    """ Imports the dbNSFP annotation table
          :param HailContext hc: The Hail context
          :param String sourcePath: Annotation table path
          :param String destinationPath: Path where the loaded annotation table will be put
          :param String nPartitions: Number of partitions
    """
    print("Annotation dbNSFP table path is " + sourcePath)
    #Variant(`#chr`,`pos(1-coor)`.toInt,`ref`,`alt`)
    table = hl.import_table(sourcePath,min_partitions=nPartitions)
    table = table.annotate(variant = hl.struct(chrom=table["#chr"],pos=table["pos(1-coor)"],ref=table["ref"],alt=table["alt"])).key_by('variant')
    # Fields renaming. Columns starting with numbers can't be selected
    table.rename({
        '1000Gp1_AF':'Gp1_AF1000',
        '1000Gp1_AC':'Gp1_AC1000',
        '1000Gp1_EUR_AF':'Gp1_EUR_AF1000',
        '1000Gp1_ASN_AF':'Gp1_ASN_AF1000',
        '1000Gp1_AFR_AF':'Gp1_AFR_AF1000',
        'ESP6500_EA_AF ':'ESP6500_EA_AF',
        'GERP++_RS':'GERP_RS'}) \
         .select(['variant',
                  'Gp1_AF1000',
                  'Gp1_EUR_AF1000',
                  'Gp1_ASN_AF1000',
                  'Gp1_AFR_AF1000',
                  'GERP_RS',
                  'MutationTaster_score',
                  'MutationTaster_pred',
                  'phyloP46way_placental',
                  'Polyphen2_HDIV_pred',
                  'Polyphen2_HVAR_score',
                  'SIFT_pred',
                  'SIFT_score',
                  'COSMIC_ID']) \
         .write(destinationPath,overwrite=True) 
    
def importDBVcf(hc, sourcePath, destinationPath, nPartitions):
    """ Imports annotations vcfs
          :param HailContext hc: The Hail context
          :param String sourcePath: Annotation vcf path
          :param String destinationPath: Path where the loaded annotation file will be put
          :param String nPartitions: Number of partitions
    """
    print("Annotation vcf source path is " + sourcePath)
    hc.import_vcf(sourcePath,min_partitions=nPartitions).write(destinationPath,overwrite=True)

def annotateVCF(hc,variants,annotationPath,destinationPath,annotations):
    """ Adds annotations to variants based on an input vds
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the annotations can be found
         :param string destinationPath: Path were the new annotated dataset can be found
         :param string annotations: Array of annotations to add to the dataset
    """
    annotationsVds = hc.read(annotationPath).split_multi()
    variants.annotate_variants_vds(annotationsVds,expr=annotations).write(destinationPath,overwrite=True)

def annotateVCFMulti(hc, variants, annotationPath, destinationPath, annotationsMulti, annotations):
    """ Adds annotations to variants that have multiallelic INFO fields.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the Clinvar annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
         :param string annotationsMulti: Array of annotations of fields that are not split when multiallelic variants are found
         :param string annotations: Array of annotations to add to the dataset
    """
    annotationsVds = hc.read(annotationPath)
    # Getting number of multiallelics
    nMultiallelics = annotationsVds.summarize().multiallelics
    annotationsVds = annotationsVds.split_multi()
    index = '0'
    # If there are multiallelics, the aIndex annotation is created by default in the dataset.
    # This is used in Hail for INFO fields which are multiallelic, since the function 'split_multi'
    # doesn't split the info field, and we need to use the aIndex in order to get the correct value.
    if nMultiallelics:
        index = 'vds.aIndex-1'
    annotationsExpr = annotationsMulti[0] % index
    for annotation in annotationsMulti[1:]:
        annotationsExpr += "," + annotation % index
    for annotation in annotations:
        annotationsExpr += "," + annotation
    variants.annotate_variants_vds(annotationsVds,expr=annotationsExpr).write(destinationPath,overwrite=True)
    
def annotateVEP(hl, variants, destinationPath, vepPath, nPartitions):
    """ Adds VEP annotations to variants.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate 
         :param string destinationPath: Path were the new annotated dataset can be found
         :param String vepPath: VEP configuration path
         :param Int nPartitions: Number of partitions 
    """
    print("Running vep")
    print("destination is "+destinationPath)
    varAnnotated = hl.vep(variants,vepPath)
    #hl.split_multi(varAnnotated) \
    varAnnotated = varAnnotated.annotate(effs=hl.map(lambda x: 
                                                     hl.struct(
                                                         gene_name=x.gene_symbol,
                                                         effect_impact=x.impact,
                                                         transcript_id=x.transcript_id,
                                                         effect=hl.str(x.consequence_terms),
                                                         gene_id=x.gene_id,
                                                         functional_class='transcript',
                                                         amino_acid_length='',
                                                         codon_change='x.hgvsc.replace(".*:","")',
                                                         amino_acid_change='x.hgvsp.replace(".*:","")',
                                                         exon_rank='x.exon',
                                                         transcript_biotype='x.biotype',
                                                         gene_coding='str(x.cds_start)'),varAnnotated.vep.transcript_consequences)) \
      .write(destinationPath,overwrite=True)

def annotateDbNSFP(hc, variants, dbnsfpPath, destinationPath):
    """ Adds dbNSFP annotations to variants.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string dbnsfpPath: Path were the dbNSFP table can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    dbnsfp = hc.read_table(dbnsfpPath)
    variants.annotate_variants_table(dbnsfp,root='va.dbnsfp') \
            .annotate_variants_expr(expr.annotationsDbNSFP()) \
            .write(destinationPath,overwrite=True)

def annotateCADD(hc, variants, annotationPath, destinationPath):
    """ Adds CADD annotations to variants.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the CADD annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    annotateVCF(hc,variants,annotationPath,destinationPath,expr.annotationsCADD())
                                   
def annotateClinvar(hc, variants, annotationPath, destinationPath):
    """ Adds Clinvar annotations to variants.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the Clinvar annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    annotateVCF(hc,variants,annotationPath,destinationPath,expr.annotationsClinvar())

def annotateDbSNP(hc, variants, annotationPath, destinationPath):
    """ Adds dbSNP annotations to variants.
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the Clinvar annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    annotateVCF(hc,variants,annotationPath,destinationPath,expr.annotationsDbSNP())
    
def annotateGnomADEx(hc, variants, annotationPath, destinationPath):
    """ Adds gnomAD Ex annotations to a dataset. 
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the GnomAD Ex annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    annotationsMulti = expr.annotationsGnomADMulti()
    annotations = expr.annotationsGnomAD()
    annotateVCFMulti(hc,variants,annotationPath,destinationPath,annotationsMulti,annotations)

def annotateExAC(hc, variants, annotationPath, destinationPath):
    """ Adds ExAC annotations to a dataset. 
         :param HailContext hc: The Hail context
         :param VariantDataset variants: The variants to annotate
         :param string annotationPath: Path were the ExAC annotation vcf can be found
         :param string destinationPath: Path were the new annotated dataset can be found
    """
    annotateVCFMulti(hc,variants,annotationPath,destinationPath,expr.annotationsExACMulti(),[])

def annotateSomatic(hl,dataset):
    annotated = dataset.transmute_entries(sample=hl.set([hl.struct(sample=dataset.s,ad=dataset.AD,dp=dataset.DP,dpstd=dataset.DPstd,gt=dataset.GT,nprogs=dataset.info.NPROGS,progs=dataset.info.PROGS)])) \
                       .drop('rsid','qual','filters','info','VAF','old_locus','old_alleles')
    return annotated.annotate_rows(samples=hl.agg.collect_as_set(annotated.sample))
