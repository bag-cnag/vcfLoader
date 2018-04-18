## Imports

from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext
from rdconnect import config, loadVCF, annotations, index, transform
from pyspark.sql.functions import lit, concat, col
from pyspark.sql.types import StringType
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
        fileName = "variantsRaw"+chrom+".vds"
        number_partitions = configuration["number_of_partitions"]
        somatic_path = utils.buildFileName(configuration["somatic_path"],chrom)
        
        if (configuration["steps"]["loadGermline"]):
            print ("step loadGermlines")
            sourceFileNameGermline = utils.buildFileName(configuration["germline_path"],chrom)
            print("sourceFilenameGermline is " + sourceFileNameGermline)
            loadVCF.importVCF(hc,sourceFileNameGermline,destination+"/loadedGermlines/"+fileName,number_partitions)
            
        if (configuration["steps"]["loadSomatic"]):
            print ("step loadSomatics")
            somatic_paths = config.readFilesList(configuration["somatic_paths"])
            dataset_paths = configuration["dataset_paths"]
            if(dataset_paths):
                dataset_paths = config.readFilesList(dataset_paths)
            loadVCF.importVCFs(hc,somatic_paths,dataset_paths,somatic_path,number_partitions)
            
        if (configuration["steps"]["annotationVEP"]):
            print ("step loadVCF")
            print ("source file is "+destination+"/loaded/"+fileName)
            annotations.annotationsVEP(hc,str(somatic_path),str(destination+"/annotatedVEP/"+fileName),configuration["vep"],number_partitions)
            #variants= hc.sqlContext.read.load("Users/dpiscia/RD-repositories/data/output/1.1.0/dataframe/chrom1")
            #annotations.VEP2(hc,variants)
        if (configuration["steps"]["loaddbNSFP"]):
            print ("step loaddbNSFP")
            annotations.importDBTable(hc,utils.buildFileName(configuration["dbNSFP_Raw"],chrom),utils.buildFileName(configuration["dnNSFP_path"],chrom),number_partitions)


        if (configuration["steps"]["loadcadd"]):
            print ("step loaddbNSFP")
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


        if (configuration["steps"]["annotatedbNSFP"]):
            print("step annotatedbNSFP")
            variants= hc.read(destination+"/annotatedVEP/"+fileName)
            annotations.annotatedbnsfp(hc,variants,utils.buildFileName(configuration["dnNSFP_path"],chrom),destination+"/annotatedVEPdbnSFP/"+fileName)

        if (configuration["steps"]["annotatecadd"]):
            print("step annotatedbcadd")
            variants= hc.read(destination+"/annotatedVEPdbnSFP/"+fileName)
            annotations.annotateVCF(hc,variants,utils.buildFileName(configuration["cadd_path"],chrom),destination+"/annotatedVEPdbnSFPCadd/"+fileName,'va.cadd = vds.info.CADD13_PHRED')

        if (configuration["steps"]["annotateclinvar"]):
            print("step annotated clinvar")
            variants = hc.read(destination+"/annotatedVEPdbnSFPCadd/"+fileName)
            annotations.annotateClinvar(hc,variants,utils.buildFileName(configuration["clinvar_path"],""),destination+"/annotatedVEPdbnSFPCaddClinvar/"+fileName)

        if (configuration["steps"]["annotateExomesGnomad"]):
            print("step annotated exomes gnomad")
            variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvar/"+fileName)
            annotations.annotateVCF(hc,variants,utils.buildFileName(configuration["exomesGnomad_path"],chrom),destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad/"+fileName,'va.gnomAD_Ex_AC =vds.info.gnomAD_Ex_AC, va.gnomAD_Ex_AF =vds.info.gnomAD_Ex_AF')

        if (configuration["steps"]["annotateWGGnomad"]):
            print("step annotated WG gnomad")
            variants= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad/"+fileName)
            annotations.annotateVCF(hc,variants,utils.buildFileName(configuration["genomesGnomad_path"],chrom),destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomad/"+fileName,'va.gnomAD_WG_AC =vds.info.gnomAD_WG_AC, va.gnomAD_WG_AF =vds.info.gnomAD_WG_AF')

        if (configuration["steps"]["transform"]):
            print ("step transform")
            # add filter ad>0 before gt collect maybe?
            grouped= hc.read(destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadWGGnomad/"+fileName)
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
            es_conf = {
                "es.mapping.id": "id",
                "es.mapping.exclude": "id",
                "es.write.operation": "upsert"
            }
            variants = sqlContext.read.load(destination+"/variants/chrom="+chrom).select("`va.predictions`","`va.populations`","`va.clinvar_filter`","`va.indel`","`va.alt`","`v.ref`","`va.pos`","`va.samples`","`va.effs`")
            variantsRN=variants.withColumnRenamed("va.predictions","predictions") \
                .withColumnRenamed("va.populations","populations") \
                .withColumnRenamed("va.indel","indel") \
                .withColumnRenamed("va.alt","alt") \
                .withColumnRenamed("v.ref","ref") \
                .withColumnRenamed("va.pos","pos") \
                .withColumnRenamed("va.samples","samples") \
                .withColumnRenamed("va.effs","effs") \
                .withColumnRenamed("va.clinvar_filter","clinvar_filter") \
                .withColumn("chrom",lit(chrom))
            id_column = concat(col("chrom").cast(StringType()), col("pos").cast(StringType()), col("ref"), col("alt"))
            variantsRN = variantsRN.withColumn("id",id_column)
            variantsRN.printSchema()
            variantsRN.write.format("org.elasticsearch.spark.sql").options(**es_conf).option("es.nodes",configuration["elasticsearch"]["host"]).option("es.port",configuration["elasticsearch"]["port"] ).save(configuration["elasticsearch"]["index_name"]+"/"+configuration["version"],mode='append')


if __name__ == "__main__":
    # Configure OPTIONS
    conf = SparkConf().setAppName(APP_NAME)
    #in cluster this will be like
    hc = hail.HailContext()
    sqlContext = SQLContext(hc.sc)
    # Execute Main functionality
    main(hc,sqlContext)
