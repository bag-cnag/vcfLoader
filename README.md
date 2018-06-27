# vcfLoader

A tool to load and annotate vcf files to ElasticSearch through Spark and Hail (https://github.com/hail-is/hail).

The annotations used and their corresponding versions are:

* VEP (v88)
* dbNSP (2.9.3) 
* CADD (1.3)
* Clinvar (20180429)
* gnomAD exome (2.0.1)
* dbSNP (build 150)
* ExAC (0.3)

In order to upload and annotate data with *vcfLoader*, the steps to follow are:

1. Install Spark 2.1.0 (currently Hail is only compatible with either Spark 2.0.2 or 2.1.0):
Download the spark distribution and add the following lines to your *.bashrc* file (with the corresponding path where Spark is downloaded):

```
export SPARK_HOME=path/to/spark
export PATH=$PATH:path/to/spark/bin
```

Run source *.bashrc* in your command line to save the changes

2. Install Hail. 

```
sudo apt-get install g++ cmake
git clone --branch 0.1 https://github.com/broadinstitute/hail.git
cd hail
./gradlew -Dspark.version=2.1.0 shadowJar
```

Add the following lines to your .bashrc file:

```
export HAIL_HOME=/path/to/hail
export PYTHONPATH="$PYTHONPATH:$HAIL_HOME/python:$SPARK_HOME/python:`echo $SPARK_HOME/python/lib/py4j*-src.zip`"
export SPARK_CLASSPATH=$HAIL_HOME/build/libs/hail-all-spark.jar
```

For more detailed instructions, go to:

https://hail.is/docs/stable/getting_started.html#building-hail-from-source

At this point, run a pyspark shell by typing 

```
pyspark --conf='spark.sql.files.openCostInBytes=53687091200' --conf='spark.sql.files.maxPartitionBytes=53687091200'
```

in your command line. The configuration parameters set are the ones needed by Hail. Once in the shell, type 

```
import hail 
```
    
to check whether the installation was done correctly.

3. Install ElasticSearch 2.4.2. You can download the zipped distribution, unzip it and run 
    
```
cd elasticsearch
./bin/elasticsearch
```
    
 Check the address *localhost:9200* in your browser to verify that is running. 

4. Install VEP. It is one of the external annotation databases used by *vcfLoader*, and needs to be installed locally:

https://www.ensembl.org/info/docs/tools/vep/index.html

Run the following commands to get the specific version used by *vcfLoader*:

```
git clone -b release/88 https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl
```  

When asked for the downloading of cache files, choose the id corresponding to *homo_sapiens hg37*.

5. Copy the configuration files *config.json* and *vep.properties* available in *Gitea*

https://172.16.10.100/gitea/platform

and adapt them to your computer paths. 

6. The configuration file specifies the paths where the annotations can be found, and *ElasticSearch* index information. 

7. In order to run the pipeline locally, download the necessary annotations and type the following command:

```
spark-submit --master local --conf='spark.sql.files.openCostInBytes=53687091200'
--conf='spark.sql.files.maxPartitionBytes=53687091200' 
--jars /home/cluengo/hail/build/libs/hail-all-spark.jar 
--executor-memory 8G 
--py-files python/dist/vcfLoader-0.1-py2.7.egg python/rdconnect/main.py 
--chrom <chromosome> --step <pipeline_steps>
```
  
Where *chromosome* is the chromosome number (MT=23, X=24 and Y=25), and *pipeline_steps* the steps to run. The latter can be specified as a single comma-separated string with the names of the steps to run. The different options are:

```
deleteIndex, createIndex, loadVCF, loaddbNSFP, loadcadd, loadclinvar, loadExomesGnomad, 
loaddbSNP, loadExAC, annotationVEP, annotateExomesGnomad, annotatedbNSFP, annotatecadd, 
annotateclinvar, annotatedbSNP, annotateExAC,groupByGenotype,transform,toElastic
```

