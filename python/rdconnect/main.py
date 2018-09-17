## Imports

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext, SparkSession
from rdconnect import config, annotations, index, transform, utils
from pyspark.sql.functions import lit
from subprocess import call
import sys, getopt
import hail

APP_NAME = "vcfLoader"
# Usage function
def usage():
    print("main.py (-c | --chrom) <chromosome_id> (-s | --step) <pipeline_step> (-n | --nchroms) <number_chromosomes_uploaded>")

# Command line arguments parser. It extracts the chromosome and the pipeline step to run
def optionParser(argv):
    chrom = ""
    step = ""
    # The number of chromosomes uploaded is only used in the counting step (to know up to
    # which chromosome to count)
    nchroms = ""
    cores = "4"
    try:
        opts, args = getopt.getopt(argv,"c:s:n:co:",["chrom=","step=","nchroms=","cores="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-c", "--chrom"):
            chrom = arg
        elif opt in ("-s", "--step"):
            step = arg
        elif opt in ("n", "--nchroms"):
            nchroms = arg
        elif opt in ("-co", "--cores"):
            cores = arg
    return chrom, nchroms, step, cores

# Main functionality. It runs the pipeline steps
def main(hc, sqlContext, configuration, chrom, nchroms, step):
    call(["ls", "-l"])

    if (chrom == "" or step == ""):
        usage()
        sys.exit(2)
        
    configuration = config.readConfig("config.json")
    destination =  configuration["destination"] + "/" + configuration["version"]
    sourceFileName = utils.buildFileName(configuration["source_path"],chrom)
    fileName = "variantsRaw" + chrom + ".vds"
    number_partitions = configuration["number_of_partitions"]
    annotations_path = configuration["annotations_path"]

    print("sourcefilename is "+sourceFileName)

    # Pipeline steps
    if ("createIndex" in step):
        print ("step to create index")
        index.create_index(configuration["elasticsearch"]["host"],configuration["elasticsearch"]["port"],configuration["elasticsearch"]["index_name"],configuration["version"],configuration["elasticsearch"]["num_shards"],configuration["elasticsearch"]["user"],configuration["elasticsearch"]["pwd"])
        
    if ("loadVCF" in step):
        print ("step loadVCF")
        rdconnect.annotations.importVCF(hc,sourceFileName,destination+"/loaded/"+fileName,number_partitions)

    if ("loaddbNSFP" in step):
        print ("step loaddbNSFP")
        annotations.importDbNSFPTable(hc,utils.buildFileName(configuration["dbNSFP_Raw"],chrom),utils.buildFileName(configuration["dnNSFP_path"],chrom),number_partitions)

    if ("loadcadd" in step):
        print ("step loadCADD")
        annotations.importDBVcf(hc,utils.buildFileName(configuration["cadd_Raw"],chrom),utils.buildFileName(configuration["cadd_path"],chrom),number_partitions)

    if ("loadclinvar" in step):
        print ("step loadclinvar")
        annotations.importDBVcf(hc,utils.buildFileName(configuration["clinvar_Raw"],""),utils.buildFileName(configuration["clinvar_path"],""),number_partitions)

    if ("loadExomesGnomad" in step):
        print ("step load exomes gnomad")
        annotations.importDBVcf(hc,utils.buildFileName(configuration["exomesGnomad_Raw"],chrom),utils.buildFileName(configuration["exomesGnomad_path"],chrom),number_partitions)

    if ("loadExAC" in step):
        print ("step load ExAC")
        annotations.importDBVcf(hc,utils.buildFileName(configuration["ExAC_Raw"],chrom),utils.buildFileName(configuration["ExAC_path"],chrom),number_partitions)
            
    if ("annotatedbNSFP" in step):
        print("step annotate dbNSFP")
        annotations_all = hc.read(annotations_path)
        annotations_dbnsfp_table = hc.read_table(utils.buildFileName(configuration["dnNSFP_path"],chrom))
        annotations_dbnsfp = hail.VariantDataset.from_table(annotations_dbnsfp_table).split_multi()
        variants = annotations.union(hc,annotations_all,annotations_dbnsfp)
        variants = variants.annotate_variants_vds(annotations_all,"va = vds")
        annotations_path = destination+"/annotateddbNSFP/"+fileName
        annotations.annotateDbNSFP(hc,variants,utils.buildFileName(configuration["dnNSFP_path"],chrom),annotations_path)

    if ("annotatecadd" in step):
        print("step annotate dbcadd")
        cadd_path = utils.buildFileName(configuration["cadd_path"],chrom)
        variants = annotations.merge_annotations(hc,annotations_path,cadd_path)
        annotations_path = destination+"/annotatedCadd/"+fileName
        annotations.annotateCADD(hc,variants,utils.buildFileName(configuration["cadd_path"],chrom),annotations_path)

    if ("annotateclinvar" in step):
        print("step annotate clinvar")
        clinvar_path = utils.buildFileName(configuration["clinvar_path"],"")
        variants = annotations.merge_annotations(hc,annotations_path,clinvar_path)
        annotations_path = destination+"/annotatedClinvar/"+fileName
        annotations.annotateClinvar(hc,variants,utils.buildFileName(configuration["clinvar_path"],""),annotations_path)

    if ("annotateExomesGnomad" in step):
        print("step annotate exomes gnomad")
        gnomad_ex_path = utils.buildFileName(configuration["exomesGnomad_path"],chrom)
        variants = annotations.merge_annotations(hc,annotations_path,gnomad_ex_path)
        annotations_path = destination+"/annotatedExGnomAD/"+fileName
        annotations.annotateGnomADEx(hc,variants,utils.buildFileName(configuration["exomesGnomad_path"],chrom),annotations_path)

    if ("annotatedbSNP" in step):
        print("step annotate dbSNP")
        dbsnp_path = utils.buildFileName(configuration["dbSNP_path"],chrom)
        variants = annotations.merge_annotations(hc,annotations_path,dbsnp_path)
        annotations_path = destination+"/annotateddbSNP/"+fileName
        annotations.annotateDbSNP(hc,variants,utils.buildFileName(configuration["dbSNP_path"],chrom),annotations_path)
        
    if ("annotateExAC" in step):
        print("step annotate ExAC")
        exac_path = utils.buildFileName(configuration["ExAC_path"],chrom)
        variants = annotations.merge_annotations(hc,annotations_path,exac_path)
        annotations_path = destination+"/annotatedExAC/"+fileName
        annotations.annotateExAC(hc,variants,utils.buildFileName(configuration["ExAC_path"],chrom),annotations_path)

    if ("annotateVEP" in step):
        print ("step annotate VEP")
        variants = hc.read(annotations_path)
        annotations.annotateVEP(hc,variants,destination+"/annotatedVEP/"+fileName,configuration["vep"],number_partitions)

    if ("annotateVariants" in step):
        print ("step annotate variants")
        variants = hc.read(destination+"/loaded/"+fileName)
        annotations = hc.read(destination+"/annotatedVEP/"+fileName)
        variants.annotate_variants_vds(annotations,expr="va=vds").write(destination+"/annotatedVariants/"+fileName)

    # Transforming step. It sets all fields to the corresponding ElasticSearch format
    if ("transform" in step):
        print ("step transform")
        annotated = hc.read(destination+"/annotatedVariants/"+fileName)
        transform.transform(annotated,destination,chrom)

    # Uploading step. It uploads all annotated variants to ElasticSearch
    if ("toElastic" in step):
        print ("step to elastic")
        es_conf = {
            "es.net.http.auth.user": configuration["elasticsearch"]["user"],
            "es.net.http.auth.pass": configuration["elasticsearch"]["pwd"],
            "es.port": configuration["elasticsearch"]["port"]
        }
        # Getting annotated variants and adding the chromosome column
        variants = sqlContext.read.load(destination+"/variants/chrom="+chrom)\
                                  .withColumn("chrom",lit(chrom))
        variants.printSchema()
        variants.write.format("org.elasticsearch.spark.sql").options(**es_conf).save(configuration["elasticsearch"]["index_name"]+"/"+configuration["version"], mode='append')

    # Counting step to check whether the number of variants in Spark corresponds to tht number of variants that
    # have been uploaded to ElasticSearch
    if ("count" in step):
        if (nchroms == ""):
            usage()
            sys.exit(2)
        count = 0
        for chrom in range(1,int(nchroms) + 1):
            variants = sqlContext.read.load(destination+"/variants/chrom=" + str(chrom))
            count += variants.count()
        print("\nTotal number of variants: " + str(count) + "\n")

    if ("compare" in step):
        if (previous == ""):
            usage()
            sys.exit(2)
        previous = sqlContext.read.load(previous+"/variants/chrom=" + str(chrom))
        current = sqlContext.read.load(destination+"/variants/chrom=" + str(chrom))
        diff = current.subtract(previous)

if __name__ == "__main__":
    # Command line options parsing
    chrom, nchroms, step, cores = optionParser(sys.argv[1:])
    main_conf = config.readConfig("config.json")
    spark_conf = SparkConf().setAppName(APP_NAME).set('spark.executor.cores',cores)
    spark = SparkSession.builder.config(conf=spark_conf).getOrCreate()
    spark.sparkContext._jsc.hadoopConfiguration().setInt("dfs.block.size",main_conf["dfs_block_size"])
    spark.sparkContext._jsc.hadoopConfiguration().setInt("parquet.block.size",main_conf["dfs_block_size"])
    hc = hail.HailContext(spark.sparkContext)
    sqlContext = SQLContext(hc.sc)
    # Execute Main functionality
    main(hc,sqlContext,main_conf,chrom,nchroms,step)
