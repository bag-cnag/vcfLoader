
import os
import requests

import rdconnect.utils as utils
from rdconnect.classException import *

def sample_index(self = None, config = None, hl = None, log = None):
	""" This function allows to update the index information for a set of 
	experiments in data-management.

	The function requires of a GenimcData object (as self) or a configuration
	pointing to a location (aka. HDFS' folder) containing a MatrixTable (HAIL).
	The experiments in the MatrixTable will be the ones which index information 
	will be updated.
	Parameters
	
	----------
	config: ConfigFile, optional
		Configuration for this step of the pipeline. If not provided or set to
		None the configuration is looked into the GenomicData in self.
	hl: context, optional
		HAIL context. If not provided or set to None the reference to the 
		module is looked into the GenomicData in self.
	log: logger, optional
		A logger to have track of the steps used in the loading process. If not
		provided or set to None the logger is looked into the GenomicData in 
		self. If no logger is in the provided nor in the GenomicData, then no
		log is performed.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = utils.check_class_and_config(self, config, hl, log)
	self.log.info('Entering updating step "DM - INDEX"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not 'process/chrom' in config.keys():
		log.warning('Provided configuration with no chromosome attached ("process/chrom"). Chromosome 21 will be used.')
		config['process/chrom'] = 21

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "index_name" filled with "{}"'.format(config['resources/elasticsearch/index_name']))

	url = 'https://' + config['applications/datamanagement/api_exp_status'].format(config['applications/datamanagement/ip'], config['applications/datamanagement/api_group'])
	headers = { 'accept': 'application/json', 
		'Content-Type': 'application/json', 
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] }
	data = "{\"dataset\":\"" + config['resources/elasticsearch/index_name'] + "\"}"
	
	# vcf = self.hl.split_multi_hts(self.hl.import_vcf(str(source_path), array_elements_required = False, force_bgz = True, min_partitions = 2))
	# full_samples = [y.get('s') for y in vcf.col.collect()]
	mtbl = self.hl.methods.read_matrix_table(source_path)
	full_samples = [ y.get('s') for y in mtbl.col.collect() ]

	self.log.info('> Experiments in loaded VCF: {}'.format(len(full_samples)))
	self.log.debug('> First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))
	self.log.debug('> Provided update content: "{}"'.format(str(data)))
	self.log.debug('> Created query URL for data-management: {}'.format(url))

	for sam in full_samples:
		q_url = url + '&experiment=' + sam
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			q_url = url + 'forceupdate=true&experiment=' + sam
			response2 = requests.post(q_url, data = data, headers = headers, verify = False)
			if response2.status_code != 200:
				self.log.error('> Query for "{}" resulted in {}. Forced update failed with {}'.format(sam, response.status_code, response2.status_code))
			else:
				self.log.warning('> Query for "{}" resulted in {}. Forced update was successful with {}'.format(sa, response.status_code, response2.status_code))
				
				# data = "{\"rawUpload\": \"pass\",		\"rawDownload\": \"pass\", \
				# 	\"receptionPipeline\": \"pass\",	\"mapping\": \"pass\", \
				# 	\"qc\": \"pass\",					\"coverage\": \"pass\", \
				# 	\"cnv\": \"waiting\",				\"variantCalling\": \"pass\", \
				# 	\"rohs\": \"pass\",					\"genomicsdb\": \"pass\", \
				# 	\"multivcf\": \"pass\",				\"hdfs\": \"pass\", \
				# 	\"es\": \"pass\",					\"in_platform\": \"pass\", \
				# 	\"dataset\":\"" + index_name + "\" \
				# }"
				# response2 = requests.post(q_url, data = data, headers = headers, verify = False)
				# if response2.status_code != 200:
				# 	accum.append((sam, response, response2))

	return self
