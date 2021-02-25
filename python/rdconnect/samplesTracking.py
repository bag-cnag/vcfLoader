
import os
import requests

import rdconnect.utils as utils
from rdconnect.classException import *
from rdconnect.classGenome import SparseMatrix
import rdconnect.setSamplesInfo as sample


# hdfs <- abans / llista de experimetns a moyure
# genomdb <- sparse_matrix
# dm compilation <- dense
# in_platform / index <- annotated m.vcf or annotated dense matrix



def sample_index(self = None, config = None, hl = None, log = None):
	""" This function allows to update the index information for a set of 
	experiments in data-management.

	The function requires of a GenimcData object (as self) or a configuration
	pointing to a location (aka. HDFS' folder) containing a MatrixTable (HAIL).
	The experiments in the MatrixTable will be the ones which index information 
	will be updated.

	This function relies on the set_experiment from setSamplesInfo.

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
	self.log.info('Entering updating step "DM - index"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	if not isSelf:
		if not 'process/chrom' in self.config.keys() or str(self.config['process/chrom']) != '21':

			self.log.warning('Provided configuration with no chromosome attached ("process/chrom") or it was not chromosome 21. Chromosome 21 will be used.')
			self.config['process/chrom'] = '21'

			gnomeADEx_annotated_file = utils.create_chrom_filename(self.config['process/destination_file'], self.config['process/chrom'])
			gnomeADEx_annotated_path = utils.create_chrom_filename(self.config['process/destination_path'], self.config['process/chrom'])
			gnomeADEx_annotated_file = utils.destination_gnomadex(os.path.join(gnomeADEx_annotated_path, self.config['resources/elasticsearch/version']), gnomeADEx_annotated_file)

			self.log.info('Since no data was provided gnomeADEx annotated set will be loaded (set: "{}"")'.format())
			self.data = hl.read_matrix_table(gnomeADEx_annotated_file)

	full_samples = [ y.get('s') for y in self.data.col.collect() ]
	self.log.debug('Total of {0} experiments were loaded'.format(len(full_samples)))

	samples.set_experiment(log = self.log, samples = full_samples, flag = 'dataset', value = self.config['resources/elasticsearch/index_name'])

	# url = 'https://' + self.config['applications/datamanagement/api_exp_status'].format(self.config['applications/datamanagement/ip'], self.config['applications/api_group'])
	# headers = { 'accept': 'application/json', 
	# 	'Content-Type': 'application/json', 
	# 	'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
	# 	'Host': self.config['applications/datamanagement/host'] }
	# data = "{\"dataset\":\"" + self.config['resources/elasticsearch/index_name'] + "\"}"
	
	# # vcf = self.hl.split_multi_hts(self.hl.import_vcf(str(source_path), array_elements_required = False, force_bgz = True, min_partitions = 2))
	# # full_samples = [y.get('s') for y in vcf.col.collect()]
	# mtbl = self.hl.methods.read_matrix_table(source_path)
	# full_samples = [ y.get('s') for y in mtbl.col.collect() ]

	# self.log.info('> Experiments in loaded VCF: {}'.format(len(full_samples)))
	# self.log.debug('> First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))
	# self.log.debug('> Provided update content: "{}"'.format(str(data)))
	# self.log.debug('> Created query URL for data-management: {}'.format(url))

	# for sam in full_samples:
	# 	q_url = url + '&experiment=' + sam
	# 	response = requests.post(q_url, data = data, headers = headers, verify = False)
	# 	if response.status_code != 200:
	# 		q_url = url + '&forceupdate=true&experiment=' + sam
	# 		response2 = requests.post(q_url, data = data, headers = headers, verify = False)
	# 		if response2.status_code != 200:
	# 			self.log.error('> Query for "{}" resulted in {} ({}). Forced update failed with {} ({})'.format(sam, response.status_code, response.text, response2.status_code, response2.text))
	# 		else:
	# 			self.log.warning('> Query for "{}" resulted in {}. Forced update was successful with {}'.format(sam, response.status_code, response2.status_code))
				
	# 			# data = "{\"rawUpload\": \"pass\",		\"rawDownload\": \"pass\", \
	# 			# 	\"receptionPipeline\": \"pass\",	\"mapping\": \"pass\", \
	# 			# 	\"qc\": \"pass\",					\"coverage\": \"pass\", \
	# 			# 	\"cnv\": \"waiting\",				\"variantCalling\": \"pass\", \
	# 			# 	\"rohs\": \"pass\",					\"genomicsdb\": \"pass\", \
	# 			# 	\"multivcf\": \"pass\",				\"hdfs\": \"pass\", \
	# 			# 	\"es\": \"pass\",					\"in_platform\": \"pass\", \
	# 			# 	\"dataset\":\"" + index_name + "\" \
	# 			# }"
	# 			# response2 = requests.post(q_url, data = data, headers = headers, verify = False)
	# 			# if response2.status_code != 200:
	# 			# 	accum.append((sam, response, response2))

	return self


def sample_in_platform(self = None, config = None, hl = None, log = None):
	""" This function allows to update the "in_platform" information for a set 
	of  experiments in data-management.

	The function requires of a GenimcData object (as self) or a configuration
	pointing to a location (aka. HDFS' folder) containing a MatrixTable (HAIL).
	The experiments in the MatrixTable will be the ones which status information 
	will be updated.

	This function relies on the set_experiment from setSamplesInfo.

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
	self.log.info('Entering updating step "DM - in_platform"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	if not isSelf:
		if not 'process/chrom' in self.config.keys() or str(self.config['process/chrom']) != '21':

			self.log.warning('Provided configuration with no chromosome attached ("process/chrom") or it was not chromosome 21. Chromosome 21 will be used.')
			self.config['process/chrom'] = '21'

			gnomeADEx_annotated_file = utils.create_chrom_filename(self.config['process/destination_file'], self.config['process/chrom'])
			gnomeADEx_annotated_path = utils.create_chrom_filename(self.config['process/destination_path'], self.config['process/chrom'])
			gnomeADEx_annotated_file = utils.destination_gnomadex(os.path.join(gnomeADEx_annotated_path, self.config['resources/elasticsearch/version']), gnomeADEx_annotated_file)

			self.log.info('Since no data was provided gnomeADEx annotated set will be loaded (set: "{}"")'.format(gnomeADEx_annotated_file))
			self.data = hl.read_matrix_table(gnomeADEx_annotated_file)

	full_samples = [ y.get('s') for y in self.data.col.collect() ]
	self.log.debug('Total of {0} experiments were loaded'.format(len(full_samples)))

	samples.set_experiment(log = self.log, samples = full_samples, flag = 'in_platform', value = 'pass')

	# url = 'https://' + self.config['applications/datamanagement/api_exp_status'].format(self.config['applications/datamanagement/ip'], self.config['applications/api_group'])
	# headers = { 'accept': 'application/json', 
	# 	'Content-Type': 'application/json', 
	# 	'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
	# 	'Host': self.config['applications/datamanagement/host'] }
	# data = "{\"in_platform\":\"pass\"}"


	# for sam in full_samples:
	# 	q_url = url + '&experiment=' + sam
	# 	response = requests.post(q_url, data = data, headers = headers, verify = False)
	# 	if response.status_code != 200:
	# 		q_url = url + '&forceupdate=true&experiment=' + sam
	# 		response2 = requests.post(q_url, data = data, headers = headers, verify = False)
	# 		if response2.status_code != 200:
	# 			self.log.error('> Query for "{}" resulted in {} ({}). Forced update failed with {} ({})'.format(sam, response.status_code, response.text, response2.status_code, response2.text))
	# 		else:
	# 			self.log.warning('> Query for "{}" resulted in {}. Forced update was successful with {}'.format(sam, response.status_code, response2.status_code))

	return self

def sample_in_genomedb(self = None, config = None, hl = None, log = None):
	""" This function allows to update the "genomedb" information for a set 
	of experiments in data-management.

	The function requires of a SparseMatrix object (as self) or a configuration
	pointing to a location (aka. HDFS' folder) containing the MatrixTable (HAIL)
	created for a SparseMatrix. The experiments in the MatrixTable will be the 
	ones which status information will be updated.

	This function relies on the set_experiment from setSamplesInfo.

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
	self.log.info('Entering updating step "DM - genomedb"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isSelf:
		if not 'process/chrom' in self.config.keys() or str(self.config['process/chrom']) != '21':
			self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')

			self.log.warning('Provided configuration with no chromosome attached ("process/chrom") or it was not chromosome 21. Chromosome 21 will be used.')
			self.config['process/chrom'] = '21'

			sparse_matrix_path = self.config['applications/combine/sparse_matrix_path']
			if sparse_matrix_path is None:
				raise NoConfigurationException('No information on "sparse_matrix_path" was provided.')

			self.log.debug('> Argument "sparse_matrix_path" filled with "{}"'.format(sparse_matrix_path))

			path_matrix = '{0}/chrom-{1}'.format(sparse_matrix_path, chrom)
			self.log.debug('Loading sparse matrix from in {0}'.format(path_matrix))
			self.data = hl.read_matrix_table(path_matrix)

	full_samples = [ y.get('s') for y in self.data.col.collect() ]
	self.log.debug('Total of {0} experiments in sparse matrix'.format(len(full_samples)))


	url = 'https://' + self.config['applications/datamanagement/api_exp_status'].format(self.config['applications/datamanagement/ip'], self.config['applications/api_group'])
	headers = { 'accept': 'application/json', 
		'Content-Type': 'application/json', 
		'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
		'Host': self.config['applications/datamanagement/host'] }
	data = "{\"genomedb\":\"pass\"}"


	for sam in full_samples:
		q_url = url + '&experiment=' + sam
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			q_url = url + '&forceupdate=true&experiment=' + sam
			response2 = requests.post(q_url, data = data, headers = headers, verify = False)
			if response2.status_code != 200:
				self.log.error('> Query for "{}" resulted in {} ({}). Forced update failed with {} ({})'.format(sam, response.status_code, response.text, response2.status_code, response2.text))
			else:
				self.log.warning('> Query for "{}" resulted in {}. Forced update was successful with {}'.format(sam, response.status_code, response2.status_code))

	return self


def samples_in_dm(self = None, config = None, hl = None, log = None):
	""" This function allows to update the "in_platform" information for a set 
	of  experiments in data-management.

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

	self, isConfig, isHl = utils.check_class_and_config(self, config, hl, log, class_to=SparseMatrix)
	self.log.info('Entering gathering step "samples_in_dm"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not 'process/chrom' in self.config.keys() or str(self.config['process/chrom']) != '21':
		self.log.warning('Provided configuration with no chromosome attached ("process/chrom") or it was not chromosome 21. Chromosome 21 will be used.')
		self.config['process/chrom'] = '21'

	smallBatch = self.config['applications/combine/sz_small_batch']
	largeBatch = self.config['applications/combine/sz_large_batch']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
	self.log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))

	self.log.debug('> Overwriting chrom to chrom-21')
	chrom = 21

	if self is None:
		sparse_path = self.config['applications/combine/sparse_matrix_path']

		# Get version of sparse matrix
		version = path.basename(path.normpath(sparse_path))
		base = sparse_path.replace(version, '')
		self.log.debug('> Detected version of sparse matrix {}'.format(version))
		self.log.debug('> Argument "sparse_path" filled with "{}"'.format(sparse_path))

		# Load sparse matrix
		try:
			self.data = hl.read_matrix_table(_name_with_chrom(sparse_path, chrom))
			self.log.info('> Sparse matrix {}/chrom-{} was loaded'.format(version, chrom))
		except:
			self.log.error('> Sparse matrix {}/chrom-{} could not be found'.format(version, chrom))
			return 

	full_samples = [ y.get('s') for y in self.data.col.collect() ]
	self.log.debug('> Number of samples in sparse matrix: {} ({}/{})'.format(len(full_samples), full_samples[ 0 ], full_samples[ -1 ]))
	self.log.debug('> First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))

	packs, n = [], 200
	for ii in range(0, len(full_samples), n):  
		packs = ','.join(full_samples[ii:ii + n])

	self.log.debug('> Data-management will be queried {} times, each time with {} experiments'.format(len(packs), n))

	url = 'https://' + self.config['applications/datamanagement/api_sm'].format(self.config['applications/datamanagement/ip'])
	headers = { 'accept': 'application/json', 
		'Content-Type': 'application/json', 
		'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
		'Host': self.config['applications/datamanagement/host'] }

	self.log.debug('> Created query URL for data-management: {}'.format(url))

	table = {}
	for ii, samlist in enumerate(packs):
		q_url = url + '?experiment=' + samlist
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			self.log.error('> Data-management returned {} ("{}") when queried with #{} batch of experiments'.format(response.status_code, response.text, ii))
			return 
		else:
			data = json.loads(resp.content)
			table.update(data)
	print(table)
	return self