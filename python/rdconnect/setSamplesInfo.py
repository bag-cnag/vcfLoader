
import json
import requests

import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData, SparseMatrix
from rdconnect.classException import *
from rdconnect.utils import check_class_and_config

from rdconnect.classLog import VoidLog

def set_experiment(self = None, config = None, hl = None, log = VoidLog(), is_playground = False, samples = [], flag = None, value = 'pass'):
	"""
	Check guidelines on (POST :: /api/statusbyexperiment/):
		https://platform.rd-connect.eu/datamanagement/api-docs/
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "set_experiment"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	isList = True
	if len(samples) == 0:
		isList = False

	if isSelf and isList:
		self.log.error('Provided both "self" (GenomicData/SparseMatrix/DenseMatrix) and "samples" (list). Only one can be provided.')
		raise Exception('Provided both "self" (GenomicData/SparseMatrix/DenseMatrix) and "samples" (list). Only one can be provided.')

	if isSelf and not isList:
		samples = [ y.get('s') for y in self.data.col.collect() ]

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	if isSelf:
		self.log.debug('>     . Argument "self" was of class {0}'.format(str(type(self))))
	self.log.debug('> Argument "samples" was set' if isList else '> Argument "samples" was not set')
	

	url = self.config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = self.config['applications/datamanagement/api_exp_status_playground'].format(
			url, self.config['applications/api_group'])
	else:
		url = self.config['applications/datamanagement/api_exp_status'].format(
			url, self.config['applications/api_group'])
	#url = '/'.join(url.split('/')[:-1])
	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
		'Host': self.config['applications/datamanagement/host'] 
	}
	if flag is None:
		raise Exception('No "flag" was provided.')
	data = "{\"" + flag + "\": \"" + value + "\"}"

	self.log.debug('> Querying {0} experiments using url "{1}"'.format(str(len(samples)), url))
	
	for ii, sam in enumerate(samples):
		q_url = url + '&forceupdate=true&experiment=' + sam
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			q_url = url + '&forceupdate=true&experiment=' + sam
			data2 = "{\"rawUpload\": \"pass\",		\"rawDownload\": \"pass\", \
				\"receptionPipeline\": \"pass\",	\"mapping\": \"pass\", \
				\"qc\": \"waiting\",				\"coverage\": \"pass\", \
				\"cnv\": \"waiting\",				\"variantCalling\": \"pass\", \
				\"rohs\": \"pass\",					\"genomicsdb\": \"waiting\", \
				\"multivcf\": \"waiting\",			\"hdfs\": \"pass\", \
				\"es\": \"waiting\",				\"in_platform\": \"waiting\", \
				\"dataset\":\"\" \
			}"
			response2 = requests.post(q_url, data = data2, headers = headers, verify = False)
			if response2.status_code != 200:
				self.log.error('Query for "{}" resulted in {} ({}). Forced update failed with {} ({})'.format(sam, response.status_code, response.text, response2.status_code, response2.text))
			else:
				q_url = url + '&experiment=' + sam
				response = requests.post(q_url, data = data, headers = headers, verify = False)
				if response.status_code != 200:
					self.log.error('After query with "forceupdate", next query for "{}" resulted in {} ({}). Forced update resulted in {} ({})'.format(sam, response.status_code, response.text, response2.status_code, response2.text))
	return self


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


# def update_dm_by_experiment(config, log = VoidLog(), experiments = [], flag = None, value = 'pass', is_playground = False):
# 	"""Function to update the data management' status of a list of experiment.
# 	"""

# 	if config is None:
# 		raise Exception('Started "update_dm_by_experiment" and no "config" was provided.')

# 	group = config['applications/datamanagement/api_group']
# 	url = config['applications/datamanagement/ip']
# 	if not url.startswith('http://') and not url.startswith('https://'):
# 		url = 'https://{0}'.format(url)

# 	if is_playground:
# 		url = config['applications/datamanagement/api_exp_status_playground'].format(url, group)
# 	else:
# 		url = config['applications/datamanagement/api_exp_status'].format(url, group)

# 	log.info('Entering step "update_dm_by_experiment"')
# 	log.debug('> Argument "flag" filled with "{}"'.format(str(flag)))
# 	log.debug('> Argument "value" filled with "{}"'.format(str(value)))

# 	if flag is None:
# 		raise Exception('Called "update_dm_by_experiment" and no flag was provided. Use "hdfs", "genomicsdb", "es", or "in_platform".')

# 	headers = { 
# 		'accept': 'application/json', 'Content-Type': 'application/json',
# 		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
# 		'Host': config['applications/datamanagement/host'] 
# 	}
# 	data = "{\"" + flag + "\": \"" + value + "\"}"

# 	for ii, sam in enumerate(experiments):
# 		q_url = url + '&experiment=' + sam
# 		response = requests.post(q_url, data = data, headers = headers, verify = False)
# 		if response.status_code != 200:
# 			self.log.error('Query #{} for experiment {} resulted in a {} code with message "{}"'.format(str(ii), q_url, str(response.status_code), response.text))
# 			print('--> Query #{} for experiment {} resulted in a {} status with content "{}"'.format(str(ii), sam, str(response.status_code), str(response.content)))



# def update_dm(initial_vcf, data_ip, data_url, data_token, field):
# 	if not field in ("genomicsdb", "hdfs", "es", "in_platform"):
# 		raise Exception("[ERROR]: (update_dm + {}) Invalid field to be updated in data data-management.".format(field))

# 	#url = "https://platform.rd-connect.eu/datamanagement/api/statusbyexperiment/?experiment="
# 	uri = "/datamanagement/api/statusbyexperiment/?forceupdate=true&experiment="
# 	url = "https://" + data_ip + uri
# 	headers = {  'Authorization': 'Token ' + data_token, "Host": data_url }
# 	data = "{\"" + field + "\":\"pass\"}"

# 	vcf = hl.split_multi_hts(hl.import_vcf(str(initial_vcf), array_elements_required = False, force_bgz = True, min_partitions = 2))
# 	full_samples = [ y.get('s') for y in vcf.col.collect() ]

# 	print('[INFO]:   . Experiments in loaded VCF: {}'.format(len(full_samples)))
# 	print('[INFO]:   . First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))
# 	print('[INFO]:   . Provided IP for data-management: {}'.format(data_ip))
# 	print('[INFO]:   . Provided URL for data-management: {}'.format(data_url))
# 	print('[INFO]:   . Provided token for data-management: {}'.format(data_token))
# 	print('[INFO]:   . Provided update content: "{}"'.format(str(data)))
# 	print('[INFO]:   . Provided field to update in data-management: {}'.format(field))
# 	print('[INFO]:   . Created query URL for data-management: {}'.format(url))

# 	for sam in full_samples:
# 		q_url = url + sam
# 		response = requests.post(q_url, data = data, headers = headers, verify = False)
# 		if response.status_code != 200:
# 			uri = "/datamanagement/api/statusbyexperiment/?forceupdate=true&experiment="
# 			q_url = "https://" + data_ip + uri
# 			response2 = requests.post(q_url, data = data, headers = headers, verify = False)
# 			if response2.status_code != 200:
# 				raise Exception('[ERROR]   . Information for sample "{}" could not be updated. Querying for second time with extended data-body'.format(sam))
