## Imports

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
from rdconnect import config, loadVCF , annotations , index , transform
from pyspark.sql.functions import lit
import hail

from rdconnect import loadVCF,utils
## CONSTANTS
from subprocess import call
APP_NAME = "My Spark Application"

##OTHER FUNCTIONS/CLASSES

## Main functionality


def main(hc,sqlContext):
    call(["ls", "-l"])

    #hc._jvm.core.vcfToSample.hello()
    configuration = config.readConfig("config.json")
    destination =  configuration["destination"] + "/" + configuration["version"]
    for chrom in configuration["chromosome"]:
        sourceFileName = utils.buildFileName(configuration["source_path"],chrom)
        print("sourcefilename is "+sourceFileName)
        fileName = "variantsRaw"+chrom+".vds"
        number_partitions = configuration["number_of_partitions"]
        source_path = destination+"/loaded/"+fileName
        annotations_path = configuration["annotations_path"]
        
        if (configuration["steps"]["loadVCF"]):
            print ("step loadVCF")
            loadVCF.importVCF(hc,sourceFileName,source_path,number_partitions)

        if (configuration["steps"]["loaddbNSFP"]):
            print ("step loaddbNSFP")
            annotations.importDBTable(hc,utils.buildFileName(configuration["dbNSFP_Raw"],chrom),utils.buildFileName(configuration["dnNSFP_path"],chrom),number_partitions)

        if (configuration["steps"]["loadcadd"]):
            print ("step loadCADD")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["cadd_Raw"],chrom),utils.buildFileName(configuration["cadd_path"],chrom),number_partitions)

        if (configuration["steps"]["loadclinvar"]):
            print ("step loadclinvar")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["clinvar_Raw"],""),utils.buildFileName(configuration["clinvar_path"],""),number_partitions)

        if (configuration["steps"]["loadExomesGnomad"]):
            print ("step load exomes gnomad")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["exomesGnomad_Raw"],chrom),utils.buildFileName(configuration["exomesGnomad_path"],chrom),number_partitions)

        if (configuration["steps"]["loadWGGnomad"]):
            print ("step load WG gnomad")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["genomesGnomad_Raw"],chrom),utils.buildFileName(configuration["genomesGnomad_path"],chrom),number_partitions)

        if (configuration["steps"]["loaddbSNP"]):
            print ("step load dbSNP")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["dbSNP_Raw"],chrom),utils.buildFileName(configuration["dbSNP_path"],chrom),number_partitions)

        if (configuration["steps"]["loadExAC"]):
            print ("step load ExAC")
            annotations.importDBvcf(hc,utils.buildFileName(configuration["ExAC_Raw"],chrom),utils.buildFileName(configuration["ExAC_path"],chrom),number_partitions)
            
        if (configuration["steps"]["annotatedbNSFP"]):
            print("step annotate dbNSFP")
            annotations_all = hc.read(annotations_path)
            annotations_dbnsfp_table = hc.read_table(utils.buildFileName(configuration["dnNSFP_path"],chrom))
            annotations_dbnsfp = hail.VariantDataset.from_table(annotations_dbnsfp_table).split_multi()
            variants = annotations.union(hc,annotations_all,annotations_dbnsfp)
            variants = variants.annotate_variants_vds(annotations_all,"va = vds")
            annotations.annotatedbnsfp(hc,variants,annotations_dbnsfp_table,destination+"/annotateddbNSFP/"+fileName)
            annotations_path = destination+"/annotateddbNSFP/"+fileName

        if (configuration["steps"]["annotatecadd"]):
            print("step annotate dbcadd")
            cadd_path = utils.buildFileName(configuration["cadd_path"],chrom)
            variants = annotations.merge_annotations(hc,annotations_path,cadd_path)
            annotations.annotatedCadd(hc,variants,cadd_path,destination+"/annotatedCadd/"+fileName)
            annotations_path = destination+"/annotatedCadd/"+fileName

        if (configuration["steps"]["annotateclinvar"]):
            print("step annotate clinvar")
            clinvar_path = utils.buildFileName(configuration["clinvar_path"],"")
            variants = annotations.merge_annotations(hc,annotations_path,clinvar_path)
            annotations.annotateClinvar(hc,variants,clinvar_path,destination+"/annotatedClinvar/"+fileName)
            annotations_path = destination+"/annotatedClinvar/"+fileName

        if (configuration["steps"]["annotateExomesGnomad"]):
            print("step annotate exomes gnomad")
            gnomad_ex_path = utils.buildFileName(configuration["exomesGnomad_path"],chrom)
            variants = annotations.merge_annotations(hc,annotations_path,gnomad_ex_path)
            annotations.annotateGnomADEx(hc,variants,gnomad_ex_path,destination+"/annotatedExGnomAD/"+fileName)
            annotations_path = destination+"/annotatedExGnomAD/"+fileName

        if (configuration["steps"]["annotateWGGnomad"]):
            print("step annotate WG gnomad")
            gnomad_wg_path = utils.buildFileName(configuration["genomesGnomad_path"],chrom)
            variants = annotations.merge_annotations(hc,annotations_path,gnomad_wg_path)
            annotations.annotateGnomADWG(hc,variants,gnomad_wg_path,destination+"/annotatedWGGnomAD/"+fileName)
            annotations_path = destination+"/annotatedWGGnomAD/"+fileName

        if (configuration["steps"]["annotatedbSNP"]):
            print("step annotate dbSNP")
            dbsnp_path = utils.buildFileName(configuration["dbSNP_path"],chrom)
            variants = annotations.merge_annotations(hc,annotations_path,dbsnp_path)
            annotations.annotatedbSNP(hc,variants,dbsnp_path,destination+"/annotateddbSNP/"+fileName)
            annotations_path = destination+"/annotateddbSNP/"+fileName

        if (configuration["steps"]["annotateExAC"]):
            print("step annotate ExAC")
            exac_path = utils.buildFileName(configuration["ExAC_path"],chrom)
            variants = annotations.merge_annotations(hc,annotations_path,exac_path)
            annotations.annotateExAC(hc,variants,exac_path,destination+"/annotatedExAC/"+fileName)
            annotations_path = destination+"/annotatedExAC/"+fileName

        if(configuration["steps"]["annotateVCF"]):
            variants = hc.read(source_path)
            annotations_vds = hc.read(utils.buildFileName(annotations_path,chrom))
            variants.annotate_variants_vds(annotations_vds,"va=vds").write(destination+"/annotated/"+fileName,overwrite=True)
            annotations_path = destination+"/annotated/"+fileName
            
        if (configuration["steps"]["annotationVEP"]):
            print ("step annotate VEP")
            annotations_path = utils.buildFileName(annotations_path,chrom)
            print ("source file is " + annotations_path)
            annotations.annotationsVEP(hc,annotations_path,str(destination+"/annotatedVEP/"+fileName),configuration["vep"],number_partitions)

        if (configuration["steps"]["groupByGenotype"]):
            print ("step groupByGenotype")
            variants = hc.read(str(destination+"/annotatedVEP/"+fileName))
            variants.annotate_variants_expr('va.samples = gs.map(g=>  {g: g, s : s}  ).collect()').write(destination+"/annotatedSamples/"+fileName,overwrite=True)

        if (configuration["steps"]["transform"]):
            print ("step transform")
            # add filter ad>0 before gt collect maybe?
            grouped= hc.read(destination+"/annotatedSamples/"+fileName)
            grouped.variants_table().to_dataframe().printSchema()
            transform.transform(grouped,destination,chrom)
            
        if (configuration["steps"]["deleteIndex"]):
            print ("step to delete index")
            index.delete_index(configuration["elasticsearch"]["host"],configuration["elasticsearch"]["port"],configuration["elasticsearch"]["index_name"],configuration["version"])

        if (configuration["steps"]["createIndex"]):
            print ("step to create index")
            index.create_index(configuration["elasticsearch"]["host"],configuration["elasticsearch"]["port"],configuration["elasticsearch"]["index_name"],configuration["version"])

        if (configuration["steps"]["toElastic"]):
            print ("step to elastic")
            variants = sqlContext.read.load(destination+"/variants_merged/chrom="+chrom).select("`va.freqInt`","`va.predictions`","`va.populations`","`va.clinvar_filter`","`va.indel`","`va.alt`","`v.ref`","`va.pos`","`va.samples`","`va.effs`")
            variantsRN=variants.withColumnRenamed("va.predictions","predictions") \
                .withColumnRenamed("va.populations","populations") \
                .withColumnRenamed("va.indel","indel") \
                .withColumnRenamed("va.alt","alt") \
                .withColumnRenamed("v.ref","ref") \
                .withColumnRenamed("va.pos","pos") \
                .withColumnRenamed("va.freqInt","freqInt") \
                .withColumnRenamed("va.samples","samples") \
                .withColumnRenamed("va.effs","effs") \
                .withColumnRenamed("va.clinvar_filter","clinvar_filter") \
                .withColumn("chrom",lit(chrom))
            variantsRN.printSchema()
            variantsRN.write.format("org.elasticsearch.spark.sql").option("es.nodes",configuration["elasticsearch"]["host"]).option("es.port",configuration["elasticsearch"]["port"] ).save(configuration["elasticsearch"]["index_name"]+"/"+configuration["version"], mode='append')

if __name__ == "__main__":
    # Configure OPTIONS
    conf = SparkConf().setAppName(APP_NAME)
    #in cluster this will be like
    hc = hail.HailContext()
    sqlContext = SQLContext(hc.sc)
    # Execute Main functionality
    main(hc,sqlContext)
