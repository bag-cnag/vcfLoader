



def importVCF(hc,source_path, destination_path,number_partitions):
    try:
        print ("reading vcf from "+ source_path )
        vcf= hc.import_vcf(str(source_path),force_bgz=True,min_partitions=number_partitions).split_multi()
        print ("writing vds to" + destination_path)
        vcf.write(destination_path,overwrite=True)
        return True
    except ValueError:
        print (ValueError)
        return "error in importing vcf"
