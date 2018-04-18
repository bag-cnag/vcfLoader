from rdconnect import annotations

def importVCF(hc,source_path, destination_path,number_partitions):
    try:
        print ("reading vcf from "+ source_path )
        vcf= hc.import_vcf(str(source_path),force_bgz=True).split_multi()
        print ("writing vds to" + destination_path)
        vcf.repartition(number_partitions).write(destination_path,overwrite=True)
        return True
    except ValueError:
        print (ValueError)
        return "error in importing vcf"

def merge(dataset, other):
    dataset_variants = dataset.drop_samples()
    dataset_variants = dataset_variants.annotate_variants_expr("va = let c = va in drop(va,samples)")
    other_variants = other.drop_samples()
    other_variants = other_variants.annotate_variants_expr("va = let c = va in drop(va,samples)")
    variants = dataset_variants.union(other_variants).deduplicate()
    variants = variants.annotate_variants_vds(dataset,"va.samples = vds.samples")
    variants = variants.annotate_variants_vds(other,"va.samples = if(isMissing(va.samples)) vds.samples else if(isMissing(vds.samples)) va.samples else va.samples.extend(vds.samples)")
    return variants

def importVCFs(hc, file_paths, dataset_paths, destination_path, num_partitions):
    nFiles = len(file_paths)
    if(nFiles > 0) :
        try:
            merged = hc.import_vcf(file_paths[0],force_bgz=True,min_partitions=num_partitions).split_multi()
            merged = annotations.annotateSomaticSamples(merged) 
            for file_path in file_paths[1:]:
                dataset = hc.import_vcf(file_path,force_bgz=True,min_partitions=num_partitions).split_multi()
                dataset = annotations.annotateSomaticSamples(dataset)
                merged = merge(merged,dataset)
            for dataset_path in dataset_paths:
                dataset = hc.read(dataset_path)
                merged = merge(merged,dataset)
            merged.write(destination_path,overwrite=True)
        except ValueError:
            print (ValueError)
            return "error in loading vcf"
    else:
        print ("Empty file list")
                            
