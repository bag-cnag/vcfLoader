import os
import sys
import getopt
import logging

import config
import moveData as mv
import loadGenetics as load
import annotateGenetics as annotate
from classGenome import add_funcs_from_module


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
def main(sqlContext, sc, main_conf, chrom, step, somaticFlag):
	var = None
	log = create_logger('vcfLoader')
	now = datetime.now()
	log.info('Staring PIPELINE at {}/{}/{} {}:{}:{}'.format(
		now.year,now.month,now.day,now.hour,now.minute,now.second)
	)

	if chrom == '':
		stop_pipeline(log, 'No chromosome was provided')
	if syep == '':
		stop_pipeline(log, 'No pipeline was provided')

	step = step.split(',')

    
	if 'move_gvcf' in step:
		mv.move_gvcf(log, chrom, 
			main_conf['process']['moving_gvcf_list'], 
			main_conf['process']['moving_from'], 
			main_conf['process']['moving_to'], 
			main_conf['process']['moving_s1'], 
			main_conf['process']['moving_s2']
		)

	if 'load_dense_matrix' in step:
		# fileName = "variants-chrom-{0}-mtx-{1}.ht".format(str(chrom), str(ii))
		destination_path = main_conf['process']['destination_path']
		destination_path = None if destination_path == '' else os.path.join(destination_path, 'dense_matrices')
		var = load.dense_matrix(log, hl, main_conf['process']['source_path'], destination_path)

    if 'VEP' in step:
		destination_path = main_conf['process']['destination_path']
		destination_path = None if destination_path == '' else destination_path
		var = annotate.vep(log, hl, var, main_conf['process']['source_path'], destination_path, 
			main_conf['annotation']['clean']['vep'], main_conf['process']['autosave'])

	if 'dbNSFP' in step:
		destination_path = main_conf['process']['destination_path']
		destination_path = None if destination_path == '' else destination_path
		var = annotate.dbNSFP(log, hl, var, main_conf['process']['source_path'], destination_path, 
			main_conf['annotation']['clean']['dbNSFP'], main_conf['process']['autosave'])




if __name__ == "__main__":
    # Command line options parsing
    chrom, path, step, cores, somaticFlag = optionParser(sys.argv[1:])
    main_conf = config.readConfig(path)
    spark_conf = SparkConf().setAppName(APP_NAME).set('spark.executor.cores',cores)
    spark = SparkSession.builder.config(conf = spark_conf).getOrCreate()
    spark.sparkContext._jsc.hadoopConfiguration().setInt("dfs.block.size", main_conf["resources"]["dfs_block_size"])
    spark.sparkContext._jsc.hadoopConfiguration().setInt("parquet.block.size", main_conf["resources"]["dfs_block_size"])
    hl.init(spark.sparkContext, tmp_dir = "hdfs://rdhdfs1:27000/test/tmp")
    sqlContext = SQLContext(hl.spark_context())
    # Execute Main functionality
    main(sqlContext, spark.sparkContext, main_conf, chrom, step, somaticFlag)
