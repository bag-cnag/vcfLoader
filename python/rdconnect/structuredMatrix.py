import sys
import json
import requests
import rdconnect.utils as utils
import rdconnect.getSamplesInfo as get

from os import path, system
from collections import Counter
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
		destination_path = destination_hdfs
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
	
	clean_to_process = []
	for item in experiments:
		if filesystem == 'ceph':
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

	#print("clean_to_process:", clean_to_process)

	try:
		self.data = hl.read_matrix_table(_name_with_chrom(sparse_path, chrom))
		self.log.info('Sparse matrix {}/chrom-{} was loaded'.format(version, chrom))
		sm_loaded = True
	except:
		self.data = None
		self.log.info('Sparse matrix {}/chrom-{} could not be found and will be created'.format(version, chrom))
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
	self.log.info('Starting step 1 - creation of cumulative matrices of {} experiments, incrementing {} experiments at a time'.format(largeBatch, smallBatch))
	batches = _create_batches(clean_to_process, version, largeBatch, smallBatch)
	

	#print("batches:", batches)

	last = None
	for idx1, batch in enumerate(batches):
		self.log.info('Processing large batch {}/{} {}'.format(idx1, len(batches), batch[ 'version' ]))

		accum = None
		for idx2, pack in enumerate(batch[ 'content' ]):
			vsr = pack[ 'version' ]
			if idx2 == len(batch[ 'content' ]) - 1:
				vsr = batch[ 'version' ]
			small_batch_path = _name_with_chrom(path.join(base, vsr), chrom)
			self.log.info('     . Loading pack #{} of {} gVCF ({})'.format(idx2, len(pack[ 'content' ]), small_batch_path))
			for f in pack['content']:
				print(f)
			last = _load_gvcf(self.hl, pack[ 'content' ], small_batch_path, accum, chrom, config[ 'applications/combine/partitions_chromosome' ])
			accum = small_batch_path

	# Collect all the small sparse matrix and iteratively accumulate them
	revisions_to_collect = [ pack[ 'version' ] for pack in batches ]
	if sm_loaded:
		revisions_to_collect = [ version ] + revisions_to_collect

	self.log.info('Starting step 2 - merging {} cumulative matrices'.format(len(revisions_to_collect)))
	print("revisions_to_collect", revisions_to_collect)
	
	for ii in range(1, len(revisions_to_collect)):
		print(revisions_to_collect[ ii-1 ], revisions_to_collect[ ii ], chrom)
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
	sm_1 = sm_1.key_rows_by('locus')
	sm_2 = sm_2.key_rows_by('locus')
	comb = combine_gvcfs([ sm_1 ] + [ sm_2 ])
	print(type(comb))
	comb.write(smD, overwrite = True)
	print(type(comb))
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
	"""
	interval = utils.get_chrom_intervals(chrom, partitions, hl)
	vcfs = [ transformFile(mt) for mt in importFiles([ x[ 'file' ] for x in experiments ]) ]

	if previous_version_path == None:
		comb = combine_gvcfs(vcfs)
	else:
		previous = hl.read_matrix_table(previous_version_path)
		previous = previous.key_rows_by('locus')
		comb = combine_gvcfs([ previous ] + vcfs)

	comb.write(version_path, overwrite = True)
	return comb
	"""

	interval = utils.get_chrom_intervals(chrom, partitions, hl)
	exp = []
	for idx, ex in enumerate(experiments):
		print("Processing file {} ({})".format(ex[ 'file' ], idx))
		x = importFiles([ ex[ 'file' ] ])
		x = transformFile(x)
		if(idx > 0):
			exp = combine_gvcfs(exp + x)
		exp.write(version_path, overwrite = True)
	
	if previous_version_path != None:
		previous = hl.read_matrix_table(previous_version_path)
		previous = previous.key_rows_by('locus')
		comb = combine_gvcfs([ previous ] + exp)
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


def append_to_dense_matrices(self = None, config = None, hl = None, log = VoidLog(), experiments = []):
	isSelf = True
	if self is None:
		isSelf = False

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

	if not isSelf:
		if sparse_matrix_path is None:
			raise NoConfigurationException('No information on "sparse_matrix_path" was provided.')

		path_matrix = '{0}/chrom-{1}'.format(sparse_matrix_path, chrom)
		self.log.debug('> Loading sparse matrix from in {}'.format(path_matrix))
		self.data = hl.read_matrix_table(path_matrix)

	experiments_in_matrix = [ x.get( 's' ) for x in self.data.col.collect() ]
	self.log.debug('> Total of {0} experiments in sparse matrix'.format(len(experiments_in_matrix)))

	#if len(experiments) == 0:
	#	self.log.info('No experiments were provided, DM will be queried to obtain the experiments to add to dense matrices (multivcf & es: waiting)')
	#	experiments = get.experiments_in_dm_traking([ (x, '') for x in experiments_in_matrix ], self.config, self.log)
	#
	#exp_sts = _get_experiments_to_dm_(self.config, self.log)
	#exp_sts = [ x['RD_Connect_ID_Experiment'] for x in exp_sts ]
	#
	#to_add = [ x for x in experiments.keys() if x in exp_sts ]
	#to_add = [ [ x ] + experiments[ x ].split('//') for x in to_add ]
	#dm_to_create = sorted(list(set([ x[ 2 ] for x in to_add ])))

	to_add = [ (x, '', '0') for x in experiments_in_matrix ]
	#to_add = [['AS5135', '', '0'], ['AS5134', '', '0'], ['AS5133', '', '0'], ['AS5137', '', '0'], ['AS5138', '', '0']]
	dm_to_create = [ '0' ]
	print("to_add:", len(to_add), to_add)
	print("dm_to_create:", dm_to_create)
	
	try:
		for idx, dm in enumerate( dm_to_create ):
			self.log.debug( "Flatting and filtering dense matrix {} (#{})".format(idx, dm))
			sam = hl.literal([ x[ 0 ] for x in to_add if x[ 2 ] == dm ], 'array<str>')
			#sam = hl.literal(experiments, 'array<str>')
			small_matrix = self.data.filter_cols(sam.contains(self.data[ 's' ]))
			small_matrix = small_matrix.key_rows_by('locus', 'alleles')
			small_matrix = hl.experimental.sparse_split_multi(small_matrix, filter_changed_loci = True)
			small_matrix = hl.experimental.densify(small_matrix)
			small_matrix = small_matrix.filter_rows(hl.agg.any(small_matrix.GT.is_non_ref()))
			small_matrix = hl.split_multi_hts(small_matrix, keep_star = False)	
			path = '{0}/chrom-{1}-mtx-{2}'.format(dense_matrix_path, chrom, dm)
			self.log.info('Writing dense matrix {} (#{}) to disk ({})'.format(dm, idx, path))
			small_matrix.write(path, overwrite = True)
			self.log.debug("Ending writing dense matrix")
	except Exception as ex:
		raise ex

	return self



def dense_matrix_grouping(self = None, config = None, hl = None, log = VoidLog(), experiments = [], N = 1000):
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = utils.check_class_and_config(self, config, hl, log, class_to=SparseMatrix)
	self.log.info('Entering gathering step "DM - dense_matrix_grouping"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')


	smallBatch = self.config['applications/combine/sz_small_batch']
	largeBatch = self.config['applications/combine/sz_large_batch']
	sparse_path = self.config['applications/combine/sparse_matrix_path']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
	self.log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))
	self.log.debug('> Argument "experiments" filled with "{}"'.format(experiments))
	self.log.debug('> Argument "sparse_path" filled with "{}"'.format(sparse_path))

	if self is None:
		self.log.info('> Since "self" is provided "experiments" will not be used.')

		if not 'process/chrom' in self.config.keys() or str(self.config['process/chrom']) != '21':
			self.log.warning('Provided configuration with no chromosome attached ("process/chrom") or it was not chromosome 21. Chromosome 21 will be used.')
			self.config['process/chrom'] = '21'
		chrom = self.config['process/chrom']

		# Get version of sparse matrix
		version = path.basename(path.normpath(sparse_path))
		base = sparse_path.replace(version, '')
		self.log.debug('> Detected version of sparse matrix {}'.format(version))
		self.log.debug('> Argument "sparse_path" filled with "{}"'.format(sparse_path))

		# Load sparse matrix
		try:
			self.data = hl.read_matrix_table(_name_with_chrom(sparse_path, chrom))
			self.log.info('Sparse matrix {}/chrom-{} was loaded'.format(version, chrom))
		except:
			self.log.error('Sparse matrix {}/chrom-{} could not be found'.format(version, chrom))
			return 

		experiments_to_proc = [ y.get('s') for y in self.data.col.collect() ]

	else:
		self.log.info('> Since "self" is not provided "experiments" will be used.')
		self.log.debug('> Total of {0} experiments where read from file'.format(len(experiments)))
		experiments_to_proc = [ x[3] for x in experiments ]

	
	self.log.debug('> Number of samples in sparse matrix: {}'.format(len(experiments_to_proc)))
	self.log.debug('> First and last sample: {} // {}'.format(experiments_to_proc[0], experiments_to_proc[len(experiments_to_proc) - 1]))

	self.log.debug('> Query DM to gather PhenoStore ids.')
	full_experiments = get.experiment_by_group(self.config, self.log)
	#print("experiments_to_proc:", experiments_to_proc)
	#print("full_experiments\n", full_experiments)

	exp_for_ps = [ (x['RD_Connect_ID_Experiment'], x['Participant_ID']) for x in full_experiments if x['RD_Connect_ID_Experiment'] in experiments_to_proc ]
	exp_and_fam = get.experiments_and_family(exp_for_ps, self.config)
	exp_in_dm = get.experiments_in_dm_traking(exp_for_ps, self.config, self.log)
	print("exp_in_dm:", exp_in_dm)
	org_dm = _experiments_with_dm_traking_(exp_and_fam, exp_in_dm, N, self.config, self.log)
	print("org_dm:", org_dm)

	with open('spark_config/dense_matrix_assignation', 'w') as fw:
		for row in org_dm:
			fw.write('\t'.join([ str(x) for x in row ]) + '\n')

	cmd = """
		cd spark_config &&
		if [ $(git status --porcelain | wc -l) -gt 0 ]; then 
			git add dense_matrix_assignation
			git commit -m "Added 'dense_matrix_assignation'."
			git push origin $gitea_branch
		fi
	"""
	system(cmd)

	return org_dm



def _experiments_with_dm_traking_(exp_and_fam, exp_in_dm, N, config, log):
	"""
		exp_and_fam: list with (RD-Connect, PhenoStore, Family)
		exp_in_dm: list with "fam//dm_idx" were fam states for family and
			dm_index number of dense matrix
		N: number of max. experiment per dense matrix
	"""

	# Obtain experiments to be added to dense matrices
	to_be_added = [ kk for kk in exp_in_dm.keys() if exp_in_dm[ kk ] == "" ]
	to_be_added = [ [ xx for xx in exp_and_fam if xx[0] == kk ][ 0 ] for kk in to_be_added ]

	# Split the information from DM according to the experiments already in dense matrices
	exp_in_dm = [ [ kk ] + exp_in_dm[ kk ].split('//') for kk in exp_in_dm.keys() ]
	for ii in range(len(exp_in_dm)):
		if len(exp_in_dm[ii]) != 4:
			exp_in_dm[ii] = exp_in_dm[ii] + ['', '']
	log.debug('nexp_in_dm: {}'.format(exp_in_dm))

	# Get which families have already members in dense matrices
	freq_fam = dict(Counter([ x[ 1 ] for x in exp_in_dm if x[ 1 ] != '' ]))
	log.debug('freq_fam: {}'.format(freq_fam))

	# Identify which experiments have other family members in dense matrices
	to_be_added = [ kk + [ kk[ 2 ] in freq_fam.keys() ] for kk in to_be_added ]
	log.debug('to_be_added: {}'.format(to_be_added))
	# If there are experiments with members in previous dense matrices get
	# those experiments so they can be moved to new dense matrices all together
	if sum([ xx[ 3 ] for xx in to_be_added]) != 0:
		log.info('There are experiments that need to be re-assigned')
		fam_dis = list(set([ xx[ 2 ] for xx in to_be_added if xx[ 3 ]]))
		log.debug('fam_dis: {}'.format(fam_dis))
		exp_dis = [ xx for xx in exp_in_dm if xx[ 1 ] in fam_dis ]
		log.debug('exp_dis: {}'.format(exp_dis))
		to_be_added += [ [ xx[ 0 ], '', xx[ 1 ] ] for xx in exp_dis ]
		to_be_added = sorted([ [ xx[ 0 ], xx[ 2 ], -1, True ] for xx in to_be_added ], key = lambda x: x[ 1 ])
	else:
		log.info('There are NO experiments that need to be re-assigned')
		to_be_added = sorted([ [ xx[ 0 ], xx[ 2 ], -1, True ] for xx in to_be_added ], key = lambda x: x[ 1 ])

	# Let's get the last dense matrix created and check if experiments can be added
	last_dm = dict(Counter([ x[ 2 ] for x in exp_in_dm if x[ 2 ] != '' ])).items()
	if len(last_dm) > 0:
		last_dm = sorted(last_dm, key = lambda x: -int(x[ 0 ]))[ 0 ]
		if last_dm[ 1 ] >= N:
			last_dm = [int(last_dm[ 0 ]) + 1, 0]
		else:
			last_dm = [int(last_dm[ 0 ]), last_dm[ 1 ]]
	else:
		last_dm = [0, 0]
	log.debug('last_dm: {}'.format(last_dm))

	# Assign all members of the same family to the last dense matrix. If the 
	# number of experiments exceeds, create a new dense matrix.
	for fam in sorted(list(set([ xx[ 1 ] for xx in to_be_added]))):
		for ii in range(len(to_be_added)):
			if to_be_added[ ii ][ 1 ] == fam:
				to_be_added[ ii ][ 2 ] = last_dm[ 0 ]
				#to_be_added[ ii ][ 3 ] = config['resources/elasticsearch/index_name'].replace('nmtrx', str(last_dm[ 0 ]))
				last_dm[ 1 ] += 1
		if last_dm[ 1 ] >= N:
			last_dm = [int(last_dm[ 0 ]) + 1, 0]

	log.debug('to_be_added: {}'.format(to_be_added))

	# Gather the information for those experiments that have not to be updated
	not_updated = [ [ xx[ 0 ], xx[ 1 ], int(xx[ 2 ]), False ] for xx in exp_in_dm if xx[ 0 ] not in [ yy[ 0 ] for yy in to_be_added ] ]
	log.debug('not_updated: {}'.format(not_updated))

	assignation = to_be_added + not_updated
	return assignation


def _get_experiments_to_dm_(config, log):
	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	url = config['applications/datamanagement/api_exp_status_list'].format(url)

	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	
	data = { "page": 1, "pageSize": 5000, 
		"fields": [ "RD_Connect_ID_Experiment", "mapping", "variantCalling", "genomicsdb", "hdfs", "es", "in_platform" ],
		"sorted": [ { "id": "RD_Connect_ID_Experiment", "desc": False} ],
		"filtered": [
			{ "id": "variantCalling", "value": "pass" }, 
			{ "id": "hdfs",           "value": "pass" },
			{ "id": "genomicsdb",     "value": "pass" },
			{ "id": "multivcf",       "value": "waiting" },
			{ "id": "es",             "value": "waiting" }
		]
	}
	log.debug('> Querying DM using URL "{0}"'.format(url))

	response = requests.post(url, json = data, headers = headers, verify = False)
	if response.status_code != 200:
		log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
		sys.exit(2)
	
	rst = response.json()
	page = rst['_meta']['total_pages']
	rst = rst['items']

	if page > 1:
		for ii in range(2, page + 1):
			data["page"] = ii
			response = requests.post(url, json = data, headers = headers, verify = False)
			if response.status_code != 200:
				log.error('Query DM for experiment list (iteration {}) resulted in a {} message'.format(str(ii), str(response.status_code)))
				sys.exit(2)
			rst.append(response.json()['items'])

	return rst

