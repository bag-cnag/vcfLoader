## Imports

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
from rdconnect import config, loadVCF , annotations , index , transform
from pyspark.sql.functions import lit, concat, col
from pyspark.sql.types import StringType
from pyspark.sql.functions import lit
import sys, getopt
import hail

from rdconnect import loadVCF,utils
## CONSTANTS
from subprocess import call
APP_NAME = "My Spark Application"

# Usage function
def usage():
    print("main.py (-c | --chrom) <chromosome_id> (-s | --step) <pipeline_step> (-n | --nchroms) <number_chromosomes_uploaded>")

# Command line arguments parser. It extracts the chromosome and the pipeline step to run
def optionParser(argv):
    chrom = ""
    step = ""
    nchroms = ""
    try:
        opts, args = getopt.getopt(argv,"c:s:n:",["chrom=","step=","nchroms="])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-c", "--chrom"):
            chrom = arg
        elif opt in ("-s", "--step"):
            step = arg
        elif opt in ("-n", "--nchroms"):
            nchroms = arg
    return chrom, nchroms, step

# Main functionality. It runs the pipeline steps
def main(argv,hc,sqlContext):
    call(["ls", "-l"])

    # Command line options parsing
    chrom, nchroms, step = optionParser(argv)
    if (chrom == "" or step == ""):
        usage()
        sys.exit(2)
        
    configuration = config.readConfig("config.json")
    destination =  configuration["destination"] + "/" + configuration["version"]
    sourceFileName = utils.buildFileName(configuration["source_path"],chrom)
    fileName = "variantsRaw"+chrom+".vds"
    number_partitions = configuration["number_of_partitions"]

    print("sourcefilename is "+sourceFileName)

    # Pipeline steps
    if ("deleteIndex" in step):
        print ("step to delete index")
        index.delete_index(configuration["elasticsearch"]["host"],configuration["elasticsearch"]["port"],configuration["elasticsearch"]["index_name"],configuration["version"])

    if ("createIndex" in step):
        print ("step to create index")
        index.create_index(configuration["elasticsearch"]["host"],configuration["elasticsearch"]["port"],configuration["elasticsearch"]["index_name"],configuration["elasticsearch"]["num_shards"],configuration["version"])
    
    if ("loadVCF" in step):
        print ("step loadVCF")
        loadVCF.importVCF(hc,sourceFileName,destination+"/loaded/"+fileName,number_partitions)

    if ("annotationVEP" in step):
        print ("step annotate VEP")
        print ("source file is "+destination+"/loaded/"+fileName)
        annotations.annotationsVEP(hc,str(destination+"/loaded/"+fileName),str(destination+"/annotatedVEP/"+fileName),configuration["vep"],number_partitions)

    if ("loaddbNSFP" in step):
        print ("step loaddbNSFP")
        annotations.importDBTable(hc,utils.buildFileName(configuration["dbNSFP_Raw"],chrom),utils.buildFileName(configuration["dnNSFP_path"],chrom),number_partitions)

    if ("loadcadd" in step):
        print ("step loadCADD")
        annotations.importDBvcf(hc,utils.buildFileName(configuration["cadd_Raw"],chrom),utils.buildFileName(configuration["cadd_path"],chrom),number_partitions)

    if ("loadclinvar" in step):
        print ("step loadclinvar")
        annotations.importDBvcf(hc,utils.buildFileName(configuration["clinvar_Raw"],""),utils.buildFileName(configuration["clinvar_path"],""),number_partitions)

    if ("loadExomesGnomad" in step):
        print ("step load exomes gnomad")
        annotations.importDBvcf(hc,utils.buildFileName(configuration["exomesGnomad_Raw"],chrom),utils.buildFileName(configuration["exomesGnomad_path"],chrom),number_partitions)

    if ("loaddbSNP" in step):
        print ("step load dbSNP")
        annotations.importDBvcf(hc,utils.buildFileName(configuration["dbSNP_Raw"],chrom),utils.buildFileName(configuration["dbSNP_path"],chrom),number_partitions)

    if ("loadExAC" in step):
        print ("step load ExAC")
        annotations.importDBvcf(hc,utils.buildFileName(configuration["ExAC_Raw"],chrom),utils.buildFileName(configuration["ExAC_path"],chrom),number_partitions)
            
    if ("annotatedbNSFP" in step):
        print("step annotate dbNSFP")
        variants= hc.read(destination+"/annotatedVEP/"+fileName)
        annotations.annotatedbnsfp(hc,variants,utils.buildFileName(configuration["dnNSFP_path"],chrom),destination+"/annotatedVEPdbnSFP/"+fileName)

    if ("annotatecadd" in step):
        print("step annotate dbcadd")
        variants= hc.read(destination+"/annotatedVEPdbnSFP/"+fileName)
        annotations.annotateVCF(hc,variants,utils.buildFileName(configuration["cadd_path"],chrom),destination+"/annotatedVEPdbnSFPCadd/"+fileName,'va.cadd = vds.info.CADD13_PHRED')

    if ("annotateclinvar" in step):
        print("step annotate clinvar")
        variants = hc.read(destination+"/annotatedVEPdbnSFPCadd/"+fileName)
        annotations.annotateClinvar(hc,variants,utils.buildFileName(configuration["clinvar_path"],""),destination+"/annotatedVEPdbnSFPCaddClinvar/"+fileName)

    if ("annotateExomesGnomad" in step):
        print("step annotate exomes gnomad")
        variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvar/"+fileName)
        annotations.annotateGnomADEx(hc,variants,utils.buildFileName(configuration["exomesGnomad_path"],chrom),destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad/"+fileName)

    if ("annotatedbSNP" in step):
        print("step annotate dbSNP")
        variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad/"+fileName)
        annotations.annotateVCF(hc,variants,utils.buildFileName(configuration["dbSNP_path"],chrom),destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomaddbSNP/"+fileName,'va.rs = vds.rsid')

    if ("annotateExAC" in step):
        print("step annotate ExAC")
        variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomaddbSNP/"+fileName)
        annotations.annotateExAC(hc,variants,utils.buildFileName(configuration["ExAC_path"],chrom),destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomaddbSNPExAC/"+fileName)
        
    if ("groupByGenotype" in step):
        print ("step groupByGenotype")
        variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomaddbSNPExAC/"+fileName)
        variants.annotate_variants_expr('va.samples = gs.map(g=>  {g: g, s : s}  ).collect()').write(destination+"/grouped/"+fileName,overwrite=True)

    if ("transform" in step):
        print ("step transform")
        grouped= hc.read(destination+"/grouped/"+fileName)
        grouped.variants_table().to_dataframe().printSchema()
        transform.transform(grouped,destination,chrom)

    if ("toElastic" in step):
        print ("step to elastic")
        es_conf = {
            "es.mapping.id": "id"
        }
        variants = sqlContext.read.load(destination+"/variants/chrom="+chrom).select("`va.freqInt`","`va.predictions`","`va.populations`","`va.clinvar_filter`","`va.gnomad_filter`","`va.indel`","`va.alt`","`v.ref`","`va.pos`","`va.samples`","`va.effs`")
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
                           .withColumnRenamed("va.gnomad_filter","gnomad_filter") \
                           .withColumn("chrom",lit(chrom))
        id_column = concat(col("chrom").cast(StringType()), lit("-"), col("pos").cast(StringType()), lit("-"), col("ref"), lit("-"), col("alt"))
        variantsRN = variantsRN.withColumn("id",id_column) 
        variantsRN.printSchema()
        variantsRN.write.format("org.elasticsearch.spark.sql").options(**es_conf).option("es.nodes",configuration["elasticsearch"]["host"]).option("es.port",configuration["elasticsearch"]["port"] ).save(configuration["elasticsearch"]["index_name"]+"/"+configuration["version"], mode='append')

    if("count" in step):
        if (nchroms == ""):
            usage()
            sys.exit(2)
        count = 0
        for chrom in range(1,int(nchroms) + 1):
            variants = sqlContext.read.load(destination+"/variants/chrom=" + str(chrom))
            count += variants.count()
        print("\nTotal number of variants: " + str(count) + "\n")

if __name__ == "__main__":
    # Configure OPTIONS
    conf = SparkConf().setAppName(APP_NAME)
    #in cluster this will be like
    hc = hail.HailContext()
    sqlContext = SQLContext(hc.sc)
    # Execute Main functionality
    main(sys.argv[1:],hc,sqlContext)
