import sys
import json
import requests
import rdconnect.utils as utils
import rdconnect.getSamplesInfo as get
from os import path
from rdconnect.classException import *
from rdconnect.classLog import VoidLog
from rdconnect.classGenome import SparseMatrix

from hail.experimental.vcf_combiner.vcf_combiner import combine_gvcfs
from hail.experimental.vcf_combiner.vcf_combiner import transform_gvcf
from hail.experimental import sparse_split_multi



"""structuredMatrix

This module contains the functions used to create a sparse matrix and to append
experiments to an already existing sparse matrix.
"""

def append_to_sparse_matrix(self = None, config = None, hl = None, log = VoidLog(), experiments = []):
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
	"""
	self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log, class_to = SparseMatrix)
	self.log.info('Entering step "append_to_sparse_matrix"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	chrom = utils.chrom_str_to_int(str(self.config['process/chrom']))
	destination_hdfs = self.config['process/moving_to_hdfs']
	destination_ceph = self.config['process/moving_to_ceph']
	sparse_path = self.config['applications/combine/sparse_matrix_path']
	filesystem = self.config['process/filesystem']

	smallBatch = self.config['applications/combine/sz_small_batch']
	largeBatch = self.config['applications/combine/sz_large_batch']


	chrom_str = chrom
	if chrom_str == '23':
		chrom_str = 'MT'
	elif chrom_str == '24':
		chrom_str = 'X'
	elif chrom_str == '25':
		chrom_str = 'Y'

	
	self.log.debug('> Argument "chrom" filled with "{}/{}"'.format(chrom, chrom_str))
	if filesystem == 'hdfs':
		self.log.debug('> Argument "destination_hdfs" filled with "{}" it will be used'.format(destination_hdfs))
		self.log.debug('> Argument "destination_ceph" filled with "{}" it will not be used'.format(destination_ceph))
		destination_path = destination_hdfds
	else:
		self.log.debug('> Argument "destination_hdfs" filled with "{}" it will not be used'.format(destination_hdfs))
		self.log.debug('> Argument "destination_ceph" filled with "{}" it will be used'.format(destination_ceph))
		destination_path = destination_ceph
	self.log.debug('> Argument "experiments" filled with "{}"'.format(experiments))
	self.log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
	self.log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))
	self.log.debug('> Argument "sparse_path" filled with "{}"'.format(sparse_path))

	# # Get experiments to load from DM
	# url = config['applications/datamanagement/ip']
	# if not url.startswith('http://') and not url.startswith('https://'):
	# 	url = 'https://{0}'.format(url)

	# url = config['applications/datamanagement/api_exp_status_list'].format(url)

	# headers = { 
	# 	'accept': 'application/json', 'Content-Type': 'application/json',
	# 	'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
	# 	'Host': config['applications/datamanagement/host'] 
	# }
	# data = "{\"page\": 1, \"pageSize\": " + str(queryBatch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"], \"sorted\":[{\"id\":\"RD_Connect_ID_Experiment\",\"desc\":false}], \"filtered\":[{\"id\":\"variantCalling\",\"value\":\"pass\"},{\"id\":\"rohs\",\"value\":\"pass\"},{\"id\":\"in_platform\",\"value\":\"waiting\"}]}"
	# self.log.debug('> Querying DM using url "{0}"'.format(url))

	# response = requests.post(url, data = data, headers = headers, verify = False)
	# if response.status_code != 200:
	# 	self.log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
	# 	sys.exit(2)

	# to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
	# self.log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))

	# all_group = get.experiment_by_group(config, self.log, False)
	# self.log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))

	# to_process = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]

	# clean_to_process = []
	# for idx, itm in enumerate(to_process):
	# 	clean_to_process.append({
	# 		'file': source_path.replace('[owner]', itm['Owner'])\
	# 			.replace('[patient-id]', itm['RD_Connect_ID_Experiment'])\
	# 			.replace('[chromosome]', str(chrom_str)),
	# 		'id': itm['RD_Connect_ID_Experiment'],
	# 		'pid': itm['Participant_ID']
	# 	})

	#experiments = ['AS5120', 'AS5121', 'AS5122', 'AS5123', 'AS5124', 'AS5125', 'AS5126', 'AS5127', 'AS5128']
	clean_to_process = []
	for item in experiments:
		if filesystem == 'chep':
			clean_to_process.append({
				'file': 's3a://cnag/' + item[1],
				'id': item[3]
			})
		else:
			clean_to_process.append({
				'file': item[2],
				'id': item[3]
			})

	# Get version of sparse matrix
	version = path.basename(path.normpath(sparse_path))
	base = sparse_path.replace(version, '')
	self.log.debug('> Detected version of sparse matrix {}'.format(version))

	try:
		self.data = hl.read_matrix_table(_name_with_chrom(sparse_path, chrom))
		self.log.info('> Sparse matrix {}/chrom-{} was loaded'.format(version, chrom))
		sm_loaded = True
	except:
		self.data = None
		self.log.info('> Sparse matrix {}/chrom-{} could not be found and will be created'.format(version, chrom))
		sm_loaded = False

	# Check for loaded experiments
	#if sm_loaded:
	#	x = [ y.get('s') for y in self.data.col.collect() ]
	#	self.log.debug('> Loaded sparse matrix contains {} experiments'.format(len(x)))
	#	y = [ z for z in clean_to_process if z['id'] in x ]
	#	if len(y) != 0:
	#		self.log.error('> {} experiments are already loaded'.format(len(y)))
	#		clean_to_process = [ z for z in clean_to_process if z['id'] not in x ]

	# Create batches of samples to be loaded
	self.log.info('> Starting step 1 - creation of cumulative matrices of {} experiments, incrementing {} experiments at a time'.format(largeBatch, smallBatch))
	batches = _create_batches(clean_to_process, version, largeBatch, smallBatch)
		
	for idx1, batch in enumerate(batches):
		self.log.info('> Processing large batch {}/{} {}'.format(idx1, len(batches), batch[ 'version' ]))

		accum = None
		for idx2, pack in enumerate(batch[ 'content' ]):
			vsr = pack[ 'version' ]
			if idx2 == len(batch[ 'content' ]) - 1:
				vsr = batch[ 'version' ]
			small_batch_path = _name_with_chrom(path.join(base, vsr), chrom)
			self.log.info('     > Loading pack #{} of {} gVCF ({})'.format(idx2, len(pack[ 'content' ]), small_batch_path))
			for f in pack['content']:
				print(f)
			_load_gvcf(self.hl, pack[ 'content' ], small_batch_path, accum, chrom, config[ 'applications/combine/partitions_chromosome' ])
			accum = small_batch_path

	# Collect all the small sparse matrix and iteratively accumulate them
	revisions_to_collect = [ pack[ 'version' ] for pack in batches ]
	revisions_to_collect = [ (self.data, version) ] + revisions_to_collect


	self.log.info('> Starting step 2 - merging {} cumulative matrices'.format(len(revisions_to_collect)))
	last = None
	for ii in range(1, len(revisions_to_collect)):
		print(ii, revisions_to_collect[ ii ])
		last = _combine_mt(self.hl, base, revisions_to_collect[ ii-1 ], revisions_to_collect[ ii ], utils.version_bump(revisions_to_collect[ ii ], 'version'), chrom)

	self.data = last
	return self


def _name_with_chrom(base, chrom):
		return path.join(base, 'chrom-{}'.format(chrom))

def _combine_mt(hl, base, ver1, ver2, verD, chrom):
	sm1 = _name_with_chrom(path.join(base, ver1), chrom)
	sm2 = _name_with_chrom(path.join(base, ver2), chrom)
	smD = _name_with_chrom(path.join(base, verD), chrom)
	print( '[_combine_mt]: merging "{}" and "{}" and saving it to "{}"'.format(sm1, sm2, smD))
	sm_1 = hl.read_matrix_table(sm1)
	sm_2 = hl.read_matrix_table(sm2)
	comb = combine_gvcfs([ sm_1 ] + [ sm_2 ])
	comb.write(smD, overwrite = True)
	return comb


def _load_gvcf(hl, experiments, version_path, previous_version_path, chrom, partitions):
	def transformFile(mt):
		x = transform_gvcf(mt.annotate_rows(
			info = mt.info.annotate(MQ_DP = hl.null(hl.tint32), VarDP = hl.null(hl.tint32), QUALapprox = hl.null(hl.tint32))
		))
		return x
	def importFiles(files):
		x = hl.import_vcfs(
			files,
			partitions = interval[ 'interval' ], 
			reference_genome = interval[ 'reference_genome' ], 
			array_elements_required = interval[ 'array_elements_required' ]
		)
		return x

	interval = utils.get_chrom_intervals(chrom, partitions, hl)
	vcfs = [ transformFile(mt) for mt in importFiles([ x[ 'file' ] for x in experiments ]) ]

	if previous_version_path == None:
		comb = combine_gvcfs(vcfs)
	else:
		previous = hl.read_matrix_table(previous_version_path)
		comb = combine_gvcfs([ previous ] + vcfs)

	comb = comb.key_rows_by('locus', 'alleles')
	comb.write(version_path, overwrite = True)


def _create_batches(experiments, version, largeSize = 500, smallSize = 100):
	""" Function to create the batches of experiments to be loaded
	and appended into the sparse matrix.
	"""
	cnt = 0
	rst = []
	lrg = []
	sml = []

	for idx, itm in enumerate(experiments):
		if len(sml) >= smallSize:
			version = utils.version_bump(version, 'iteration')
			cnt += smallSize
			lrg.append({ 'version': version, 'content': sml })
			sml = []
		
		if cnt >= largeSize:
			version = utils.version_bump(version, 'revision')
			rst.append({ 'version': version, 'content': lrg })
			lrg = []
			cnt = 0
		
		sml.append(itm)

	if len(sml) > 0:
		version = utils.version_bump(version, 'iteration')
		lrg.append({ 'version': version, 'content': sml })
		version = utils.version_bump(version, 'revision')
		rst.append({ 'version': version, 'content': lrg })

	return rst

def append_to_dense_matrices(self = None, config = None, hl = None, log = VoidLog()):
	self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log, class_to = SparseMatrix)
	self.log.info('Entering step "append_to_dense_matrices"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')


	chrom = utils.chrom_str_to_int(str(config['process/chrom']))
	dense_matrix_path = self.config['applications/combine/dense_matrix_path']
	sparse_matrix_path = self.config['applications/combine/sparse_matrix_path']
	sz_large_batch = self.config['applications/combine/sz_large_batch']

	self.log.debug('> Argument "chrom" filled with "{}"'.format(chrom))
	self.log.debug('> Argument "dense_matrix_path" filled with "{}"'.format(dense_matrix_path))
	self.log.debug('> Argument "sparse_matrix_path" filled with "{}"'.format(sparse_matrix_path))

	#mapping = load_table_log(sq, '{0}/mapping'.format(dense_matrix_path))

	if sparse_matrix_path is None:
		raise NoConfigurationException('No information on "sparse_matrix_path" was provided.')

	path_matrix = '{0}/chrom-{1}'.format(sparse_matrix_path, chrom)
	self.log.debug('Loading sparse matrix from in {0}'.format(path_matrix))
	sparse_matrix = hl.read_matrix_table(path_matrix)

	experiments_in_matrix = [ x.get( 's' ) for x in sparse_matrix.col.collect() ]
	self.log.debug('Total of {0} experiments in sparse matrix'.format( len( experiments_in_matrix ) ))

	idx = 0
	try:
		#for idx, batch in enumerate( mapping ):
		#	self.log.debug( "Flatting and filtering dense matrix {0} (sz: {1}) --> {2} - {3}".format( idx, len( batch ), batch[0], batch[len(batch) - 1] ) )
		#	sam = hl.literal([ x[ 0 ] for x in batch ], 'array<str>')
		sam = hl.literal(experiments_in_matrix, 'array<str>')
		small_matrix = sparse_matrix.filter_cols(sam.contains(sparse_matrix[ 's' ]))
		small_matrix = hl.experimental.densify(small_matrix)
		small_matrix = small_matrix.filter_rows(hl.agg.any(small_matrix.LGT.is_non_ref()))
		small_matrix = hl.split_multi_hts(small_matrix, keep_star = False)	
		path = '{0}/chrom-{1}-mtx-{2}'.format(dense_matrix_path, chrom, idx)
		self.log.info('Writing dense matrix {} to disk ({})'.format(idx, path))
		small_matrix.write(path, overwrite = True)
		self.log.debug( "Ending writing dense matrix" )
	except Exception as ex:
		raise ex

	return self



def dense_matrix_grouping(self = None, config = None, hl = None, log = VoidLog()):
	self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log, class_to = SparseMatrix)
	self.log.info('Entering step "dense_matrix_grouping"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	self.log.debug('OVERWRITING chrom to chrom-21')

	chrom = 21
	sparse_path = path.join(self.config['applications/combine/sparse_matrix_path'], 'chrom-{}'.format(chrom))
	sparse_matrix = hl.read_matrix_table(sparse_path)

	experiments_in_matrix = [ x.get( 's' ) for x in sparse_matrix.col.collect() ]
	self.log.debug('Obtained a total of {} experiments from sparse matrix (chrom 21)'.format(len(experiments_in_matrix)))

	all_group = get.experiment_by_group(config, self.log, False)
	self.log.debug('Obtained a total of {} experiments for the group'.format(len(all_group)))

	full_ids_in_matrix = [ x for x in all_group if x[ 'RD_Connect_ID_Experiment' ] in experiments_in_matrix ]
	print('full_ids_in_matrix', len( full_ids_in_matrix ))
	print('\t', full_ids_in_matrix[ : 10 ])

	experiments_and_families = get.experiments_and_family(full_ids_in_matrix, self.config)
	print('experiments_and_families', len( experiments_and_families ))
	print('\t', experiments_and_families[ : 10 ])


	# # Get all the experiments that have to processed from data-management
	# experiments_in_group = getExperimentByGroup( group, url_project, host_project, token, prefix_hdfs, chrom, max_items_batch, is_playground )
	# print('experiments_in_group', len( experiments_in_group ))
	# print('\t', experiments_in_group[ : 2 ])
	# full_ids_in_matrix = [ x for x in experiments_in_group if x[ 'RD_Connect_ID_Experiment' ] in experiments_in_matrix ]
	# print('full_ids_in_matrix', len( full_ids_in_matrix ))
	# print('\t', full_ids_in_matrix[ : 2 ])
	# experiments_and_families = getExperimentsByFamily( full_ids_in_matrix, url_project, gpap_id, gpap_token )
	# print('experiments_and_families', len( experiments_and_families ))

	# # Relocate experiments with no family
	# none_detected = False
	# x = len( list( set( [ x[ 2 ] for x in experiments_and_families ] ) ) )
	# for ii in range( len( experiments_and_families ) ):
	#     if experiments_and_families[ ii ][ 2 ] == '---':
	#         none_detected = True
	#         experiments_and_families[ ii ][ 2 ] = experiments_and_families[ ii ][ 0 ]
	# y = len( list( set( [ x[ 2 ] for x in experiments_and_families ] ) ) )
	# if none_detected:
	#     warnings.warn( 'Provided experiment ids got no family assigned. RD-Connect ID used as family ID for those experiments. Original families were of {} while after update are of {}.'.format( x, y ) )
	# experiments_and_families.sort(key=lambda x: x[ 0 ])

	# batches = create_batches_by_family( experiments_and_families, 1000 )
	# lgr.debug( 'Created {} batches'.format( len( batches ) ) )
	# #for ii, bat in enumerate(batches):
	# #    print('\tBatch {0}: {1} --> {2} - {3}'.format( ii, len( bat ), bat[0], bat[len(bat) - 1]))
	# for sam in batches:
	#     print(sam[0], "\t", sam[1], "\t", sam[2], "\t", sam[3])

	# log_path = '{0}/mapping'.format(dense_matrix_path)
	# rdd = sc.parallelize(batches)
	# experiments = rdd.map(lambda x: Row( RD_Connect_ID = x[ 0 ], PhenoTips = x[ 1 ], Family = x[ 2 ], DMatrix = x[ 3 ]))
	# df = sq.createDataFrame(experiments)
	# df.repartition(1).write.format('csv').mode('overwrite').save(log_path, header = 'true')
