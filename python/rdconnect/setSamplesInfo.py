
import json
import requests

import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData, SparseMatrix
from rdconnect.classException import *
from rdconnect.utils import check_class_and_config

from rdconnect.classLog import VoidLog

def set_experiment(self = None, config = None, hl = None, log = VoidLog(), is_playground = False, samples = [], flag = None, value = 'pass'):
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "ClinVar"')

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
	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
		'Host': self.config['applications/datamanagement/host'] 
	}
	data = "{\"" + flag + "\": \"" + value + "\"}"

	self.log.debug('> Querying {0} experiments using url "{1}"'.format(str(len(samples)), url))
	
	for sam in samples:
		q_url = url + sam
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			print(response.status_code)
			print(response.content)
			raise Exception("no 200")

	return self


def update_dm(initial_vcf, index_name, data_ip, data_url, data_token, field):
	if not field in ("genomicsdb", "hdfs", "es", "in_platform"):
		raise Exception("[ERROR]: (update_dm + {}) Invalid field to be updated in data data-management.".format(field))

	#url = "https://platform.rd-connect.eu/datamanagement/api/statusbyexperiment/?experiment="
	uri = "/datamanagement/api/statusbyexperiment/?forceupdate=true&experiment="
	url = "https://" + data_ip + uri
	headers = {  'Authorization': 'Token ' + data_token, "Host": data_url }
	data = "{\"" + field + "\":\"pass\"}"

	vcf = hl.split_multi_hts(hl.import_vcf(str(initial_vcf), array_elements_required = False, force_bgz = True, min_partitions = 2))
	full_samples = [ y.get('s') for y in vcf.col.collect() ]

	print('[INFO]:   . Experiments in loaded VCF: {}'.format(len(full_samples)))
	print('[INFO]:   . First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))
	print('[INFO]:   . Provided IP for data-management: {}'.format(data_ip))
	print('[INFO]:   . Provided URL for data-management: {}'.format(data_url))
	print('[INFO]:   . Provided token for data-management: {}'.format(data_token))
	print('[INFO]:   . Provided update content: "{}"'.format(str(data)))
	print('[INFO]:   . Provided field to update in data-management: {}'.format(field))
	print('[INFO]:   . Created query URL for data-management: {}'.format(url))

	for sam in full_samples:
		q_url = url + sam
		response = requests.post(q_url, data = data, headers = headers, verify = False)
		if response.status_code != 200:
			uri = "/datamanagement/api/statusbyexperiment/?forceupdate=true&experiment="
			q_url = "https://" + data_ip + uri
			response2 = requests.post(q_url, data = data, headers = headers, verify = False)
			if response2.status_code != 200:
				raise Exception('[ERROR]   . Information for sample "{}" could not be updated. Querying for second time with extended data-body'.format(sam))
