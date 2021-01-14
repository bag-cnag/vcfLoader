import json
import requests
import rdconnect.utils as utils
import rdconnect.getSamplesInfo as get
from os import path
from rdconnect.classException import *
from rdconnect.classLog import VoidLog


"""structuredMatrix

This module contains the functions used to create a sparse matrix and to append
experiments to an already existing sparse matrix.
"""

def append_to_sparse_matrix(self = None, config = None, hl = None, log = VoidLog(), largeBatch = 500, smallBatch = 100):
	""" [...]


	process/moving_to
	applications/combine/sparse_matrix_path

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

	chrom = utils.chrom_str_to_int(str(config['process/chrom']))
	source_path = self.config['process/moving_to']
	sparse_path = self.config['applications/combine/sparse_matrix_path']

	chrom_str = chrom
	if chrom_str == '23':
		chrom_str = 'MT'
	elif chrom_str == '24':
		chrom_str = 'X'
	elif chrom_str == '25':
		chrom_str = 'Y'

	self.log.debug('> Argument "chrom" filled with "{}/{}"'.format(chrom, chrom_str))
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
	self.log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))
	self.log.debug('> Argument "sparse_path" filled with "{}"'.format(sparse_path))

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
	self.log.debug('> Querying DM using url "{0}"'.format(url))

	response = requests.post(url, data = data, headers = headers, verify = False)
	if response.status_code != 200:
		self.log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
		sys.exit(2)

	to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
	self.log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))

	all_group = get.experiment_by_group(config, self.log, False)
	self.log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))

	to_process = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]

	clean_to_process = []
	for idx, itm in enumerate(to_process):
		clean_to_process.append({
			'file': source_path.replace('[owner]', itm['Owner'])\
				.replace('[patient-id]', itm['RD_Connect_ID_Experiment'])\
				.replace('[chromosome]', str(chrom_str)),
			'id': itm['RD_Connect_ID_Experiment'],
			'pid': itm['Participant_ID']
		})

	# Get version of sparse matrix
	version = path.basename(path.normpath(sparse_path)).split('.')
	self.log.debug('> Detected version of sparse matrix {}.{}.{}'.format(version[0], version[1], version[2]))

	try:
		sm = hl.read_matrix_table(os.path.join(sparse_path), 'chrom-{}'.format(chrom))
		self.log.info('> Sparse matrix {}.{}.{} was loaded'.format(version[0], version[1], version[2]))
		sm_loaded = True
	except:
		self.log.info('> Sparse matrix {}.{}.{} could not be found and will be created'.format(version[0], version[1], version[2]))
		sm_loaded = False
	

	
	print('-' * 25)
	print(clean_to_process[0])
	print('-' * 25)

	batches = _create__batches(clean_to_process, largeBatch, smallBatch)

	print(len(batches))

	for idx, elm in enumerate(batches):
		print(idx)
		for idx2, elm2 in elm['large_batch']:
			print(idx2, len(elm2))

	print("=" * 25)
	print(batches)





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
			#largeBatch.append( { 'uri': uri, 'small_batch': smallBatch } )
			largeBatch.append( { 'small_batch': smallBatch } )
			cnt += smallSize
			smallBatch = []
			added = True

		if cnt >= largeSize:
			#rst.append( { 'uri': uri, 'large_batch': largeBatch } )
			rst.append( { 'large_batch': largeBatch } )
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

		smallBatch.append( itm )

	if len( smallBatch ) != 0:
		if not bumpRev:
			#uri = utils.version_bump( uri, 'revision' )
			pass
		#rst.append( { 'uri': uri, 'large_batch': [ { 'uri': uri, 'small_batch': smallBatch } ] } )
		rst.append( { 'large_batch': [ { 'small_batch': smallBatch } ] } )
	return rst