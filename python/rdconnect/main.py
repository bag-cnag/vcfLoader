import os
import sys
import getopt
import logging
import hail as hl

from datetime import datetime
from pyspark import SparkConf, SparkContext
from pyspark.sql import SQLContext, SparkSession

import rdconnect.moveData as mv
#import rdconnect.importGenetics as cd
import rdconnect.annotateGenetics as annotate
from rdconnect.classConfig import ConfigFile
from rdconnect.classGenome import add_funcs_from_module


# USAGE SECTION ---------------------------------------------------------------
def usage():
	print("main.py (-c | --chrom) <chromosome_id> (-s | --step) <pipeline_step> (-p | --path) <config_path> (-d | --somatic_data)")

def optionParser(argv):
	chrom = ""
	nchroms = ""
	step = ""
	cores = "1"
	path = "config.json"
	somaticFlag = False
	try:
		opts, args = getopt.getopt(argv, "c:p:s:n:co:d:", ["chrom=", "path=", "step=", "nchroms=", "cores=", "somatic_data="])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-c", "--chrom"):
			chrom = arg
		elif opt in ("-p", "--path"):
			path = arg
		elif opt in ("-s", "--step"):
			step = arg
		elif opt in ("-co", "--cores"):
			cores = arg
		elif opt in ("n", "--nchroms"):
			nchroms = arg
		elif opt in ("-d", "--somatic_data"):
			if arg.lower() == 'yes':
				somaticFlag = True
			else:
				somaticFlag = False
	return chrom, path, step, cores, somaticFlag


# CONFIGURATION SECTION -------------------------------------------------------
APP_NAME = "vcfLoader"

def create_logger(name, path = None):
	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	logger = logging.getLogger(name)
	logger.setLevel(logging.DEBUG)
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	ch.setFormatter(formatter)
	logger.addHandler(ch)
	if path is not None:
		date_time = datetime.now().strftime("%y%m%d_%H%M%S")
		fh = logging.FileHandler(os.path.join(path, 'vcfLoader_{}_debug.log'.format(date_time)))
		fh.setLevel(logging.DEBUG)
		fh.setFormatter(formatter)
		logger.addHandler(fh)
	return logger


def stop_pipeline(log, msg):
	log.error('An error happened: "{}"'.format(msg))
	now = datetime.now()
	log.info('Stopping the PIPELINE at {}/{}/{} {}:{}:{}'.format(
		now.year,now.month,now.day,now.hour,now.minute,now.second)
	)
	sys.exit(2)

# MAIN METHOD SECTION ---------------------------------------------------------
def main(sqlContext, sc, config, chrom, step, somaticFlag):
	var = None
	log = create_logger('vcfLoader')
	now = datetime.now()
	log.info('Staring PIPELINE at {}/{}/{} {}:{}:{}'.format(
		now.year,now.month,now.day,now.hour,now.minute,now.second)
	)

	if chrom == '':
		stop_pipeline(log, 'No chromosome was provided')
	if step == '':
		stop_pipeline(log, 'No pipeline was provided')

	step = step.split(',')
	config.overwrite('process/chrom', chrom)
	destination_path = config['process/destination_path']

	if 'move_gvcf' in step:
		mv.gvcf(config, log)

	if 'load_dense_matrix' in step:
		#fileName = "variants-chrom-{0}-mtx-{1}.ht".format(str(chrom), str(ii))
		#destination_path = None if destination_path == '' else os.path.join(destination_path, 'dense_matrices')
		#config.overwrite('process/filename', )
		#var = load.dense_matrix(local_conf log, hl, config['process']['source_path'], destination_path)
		pass

	if 'annotateVEP' in step:
		var = annotate.vep(None, config, hl, log)

	if 'annotatedbNSFP' in step:
		var = annotate.dbNSFP(None, config, hl, log)

	if 'annotateFullDenseMatrix' in step:
		add_funcs_from_module(annotate)
		local = config.overwrite('process/autosave', False)
		var = annotate.vep(None, local, hl, log)\
				.annotate.dbNSFP(config, hl, log)




if __name__ == "__main__":
	# Command line options parsing
	chrom, path, step, cores, somaticFlag = optionParser(sys.argv[1:])
	config = ConfigFile(path)
	spark_conf = SparkConf().setAppName(APP_NAME).set('spark.executor.cores',cores)
	spark = SparkSession.builder.config(conf = spark_conf).getOrCreate()
	spark.sparkContext._jsc.hadoopConfiguration().setInt("dfs.block.size", config["resources/dfs_block_size"])
	spark.sparkContext._jsc.hadoopConfiguration().setInt("parquet.block.size", config["resources/dfs_block_size"])
	hl.init(spark.sparkContext, tmp_dir = "hdfs://rdhdfs1:27000/test/tmp")
	sqlContext = SQLContext(hl.spark_context())
	# Execute Main functionality
	main(sqlContext, spark.sparkContext, config, chrom, step, somaticFlag)
