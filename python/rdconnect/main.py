## Imports

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext, SparkSession
from rdconnect import config, annotations, index, transform, utils
from pyspark.sql.functions import lit
from subprocess import call
import sys, getopt
import hail as hl

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
def main(hl, sqlContext, configuration, chrom, nchroms, step):
    call(["ls", "-l"])

    if (chrom == "" or step == ""):
        usage()
        sys.exit(2)
        
    configuration = config.readConfig("config.json")
    destination =  configuration["destination"] + "/" + configuration["version"]
    sourceFileName = utils.buildFileName(configuration["source_path"],chrom)
    fileName = "variantsRaw" + chrom + ".ht"
    number_partitions = configuration["number_of_partitions"]

    print("sourcefilename is "+sourceFileName)

    if ("loadSomatic" in step):
        print ("step loadSomatics")
        # Read somatic vcf file
        somatic_paths = config.readFilesList(configuration["somatic_paths"])
        loaded_path = destination+"/loadedSomatic/"+fileName
        # Import and merge somatic files
        annotations.importSomatic(hl,somatic_paths,destination+"/loadedSomatic/"+fileName,number_partitions)

    if ("loaddbNSFP" in step):
        print ("step loaddbNSFP")
        annotations.importDbNSFPTable(hl,utils.buildFileName(configuration["dbNSFP_Raw"],chrom),utils.buildFileName(configuration["dnNSFP_path"],chrom),number_partitions)
        
    if ("annotateVEP" in step):
        print("Step annotate VEP")
        variants = hl.read_table(destination+"/loadedSomatic/"+fileName)
        annotations.annotateVEP(hl,variants,destination+"/annotatedVEP/"+fileName,configuration["vep"],number_partitions)

    if ("annotatedbNSFP" in step):
        print("Step annotate dbNSFP")
        variants = hc.read(destination+"/annotatedVEP/"+fileName)
        annotations.annotateDbNSFP(hl,variants,utils.buildFileName(configuration["dnNSFP_path"],chrom),destination+"/annotatedVEPdbnSFP/"+fileName)

if __name__ == "__main__":
    # Command line options parsing
    chrom, nchroms, step, cores = optionParser(sys.argv[1:])
    main_conf = config.readConfig("config.json")
    spark_conf = SparkConf().setAppName(APP_NAME).set('spark.executor.cores',cores)
    spark = SparkSession.builder.config(conf=spark_conf).getOrCreate()
    spark.sparkContext._jsc.hadoopConfiguration().setInt("dfs.block.size",main_conf["dfs_block_size"])
    spark.sparkContext._jsc.hadoopConfiguration().setInt("parquet.block.size",main_conf["dfs_block_size"])
    hl.init(spark.sparkContext)
    sqlContext = SQLContext(hl.spark_context())
    # Execute Main functionality
    main(hl,sqlContext,main_conf,chrom,nchroms,step)
