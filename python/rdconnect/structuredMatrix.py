import json
import requests
import rdconnect.utils as utils
import rdconnect.getSamplesInfo as get
from rdconnect.classException import *
from rdconnect.classLog import VoidLog


"""structuredMatrix

This module contains the functions used to create a sparse matrix and to append
experiments to an already existing sparse matrix.
"""

def append_to_sparse_matrix(self = None, config = None, hl = None, log = VoidLog(), largeBatch = 500, smallBatch = 100):
	""" [...]


	process/moving_to
	combine/sparse_matrix_path

	'applications/datamanagement/ip'
	'applications/datamanagement/api_exp_status_list'

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
	largeBatch
	smallBatch
	"""
	self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log)
	self.log.info('Entering step "append_to_sparse_matrix"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	chrom = config['process/chrom']
	chrom = utils.chrom_str_to_int(str(chrom))
	source_path = self.config['process/moving_to']
	sparse_path = self.config['combine/sparse_matrix_path']

	log.debug('> Argument "chrom" filled with "{}"'.format(chrom))
	log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
	log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))

	# Get experiments to load from DM

	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	url = config['applications/datamanagement/api_exp_status_list'].format(url)

	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	data = "{\"page\": 1, \"pageSize\": " + str(largeBatch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"], \"sorted\":[{\"id\":\"RD_Connect_ID_Experiment\",\"desc\":false}], \"filtered\":[{\"id\":\"variantCalling\",\"value\":\"pass\"},{\"id\":\"rohs\",\"value\":\"pass\"},{\"id\":\"in_platform\",\"value\":\"waiting\"}]}"
	log.debug('> Querying DM using url "{0}"'.format(url))

	response = requests.post(url, data = data, headers = headers, verify = False)
	if response.status_code != 200:
		log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
		sys.exit(2)

	to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
	log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))

	all_group = get.experiment_by_group(config, log, False)
	log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))

	to_process_group = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]

	print("hello! I got {} experiments to process from the total of {} I asked".format(len(to_process_group), largeBatch))


	print('-' * 25) 
	print(to_process[0])
	print('-' * 25)
	print(all_group[0])
	print('-' * 25)
	print(to_process_group[0])
	print('-' * 25)

	batches = _create__batches(to_process_group, largeBatch, smallBatch)

	print(len(batches))
	print(len(batches[0]))





def _create__batches(list_experiments, largeSize = 500, smallSize = 100):
	""" Function to create the batches of experiments to be loaded
	and appended into the sparse matrix.
	"""
	cnt = 0
	rst = []

	smallBatch = []
	largeBatch = []
	added = False
	bumpRev = False

	for idx, itm in enumerate( list_experiments ):   
		if len( smallBatch ) >= smallSize:
			#largeBatch.append( { 'uri': uri, 'batch': smallBatch } )
			largeBatch.append( { 'batch': smallBatch } )
			cnt += smallSize
			smallBatch = []
			added = True

		if cnt >= largeSize:
			#rst.append( { 'uri': uri, 'batches': largeBatch } )
			rst.append( { 'batches': largeBatch } )
			largeBatch = [ ]
			cnt = 0

		if added:
			if cnt + smallSize >= largeSize:
				#uri = utils.version_bump( uri, 'revision' )
				bumpRev = True
			else:
				#uri = utils.version_bump( uri, 'iteration' )
				bumpRev = False
			added = False

		smallBatch.append( item )

	if len( smallBatch ) != 0:
		if not bumpRev:
			#uri = utils.version_bump( uri, 'revision' )
			pass
		#rst.append( { 'uri': uri, 'batches': [ { 'uri': uri, 'batch': smallBatch } ] } )
		rst.append( { 'batches': [ { 'batch': smallBatch } ] } )
	return rst