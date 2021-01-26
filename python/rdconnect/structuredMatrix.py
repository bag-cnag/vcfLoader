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



"""structuredMatrix

This module contains the functions used to create a sparse matrix and to append
experiments to an already existing sparse matrix.
"""

def append_to_sparse_matrix(self = None, config = None, hl = None, log = VoidLog(), queryBatch = 500, largeBatch = 500, smallBatch = 100):
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
	self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log, class_to = SparseMatrix)
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
	self.log.debug('> Argument "queryBatch" filled with "{}"'.format(queryBatch))
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

	to_process = [ 
		"EPR002042", "EPR003231", "EPR018136", "EPR029251", "EPR029678", 
		"EPR044671", "EPR049573", "EPR051965", "EPR063942", "EPR079103", 
		"EPR084663", "EPR107689", "EPR113523", "EPR121024", "EPR153765", 
		"EPR167258", "EPR181750", "EPR182917", "EPR183701", "EPR184661", 
		"EPR185405", "EPR193852", "EPR216658", "EPR223743", "EPR263520", 
		"EPR285612", "EPR291990", "EPR305671", "EPR325618", "EPR360408", 
		"EPR370044", "EPR388888", "EPR400684", "EPR422220", "EPR422766", 
		"EPR430999", "EPR441186", "EPR441415", "EPR445201", "EPR457507", 
		"EPR462646", "EPR462953", "EPR465551", "EPR465945", "EPR468980", 
		"EPR483682", "EPR490261", "EPR491136", "EPR491842", "EPR495233", 
		"EPR503930", "EPR515435", "EPR517291", "EPR519505", "EPR520853", 
		"EPR526553", "EPR532501", "EPR543983", "EPR550773", "EPR555002", 
		"EPR567596", "EPR585246", "EPR604034", "EPR623120", "EPR662811", 
		"EPR683873", "EPR684728", "EPR688979", "EPR693989", "EPR695731", 
		"EPR702761", "EPR713819", "EPR727555", "EPR729299", "EPR754736", 
		"EPR756268", "EPR770222", "EPR770868", "EPR791504", "EPR799473", 
		"EPR815683", "EPR843651", "EPR846508", "EPR858967", "EPR870253", 
		"EPR872915", "EPR884145", "EPR889407", "EPR889571", "EPR900812", 
		"EPR905200", "EPR907717", "EPR934213", "EPR949683", "EPR988173", 
		"EPR988650", "EPR023120", "EPR167786", "EPR248369", "EPR342942", 
		"EPR363044", "EPR494310", "EPR802952", "EPR825915", "EPR857658", 
		"EPR895006", "EPR940330", "EPR970855", "EPR021681", "EPR148347", 
		"EPR356647", "EPR445362", "EPR681606", "EPR032694", "EPR033734", 
		"EPR068469", "EPR106784", "EPR120649", "EPR166943", "EPR435419", 
		"EPR456985", "EPR462421", "EPR506987", "EPR523069", "EPR559975", 
		"EPR641688", "EPR707076", "EPR748165", "EPR770552", "EPR889610", 
		"EPR935190", "EPR964066", "EPR995111", "EPR085535", "EPR090109", 
		"EPR136882", "EPR154786", "EPR234746", "EPR260721", "EPR280622", 
		"EPR315325", "EPR480495", "EPR486307", "EPR489688", "EPR546495", 
		"EPR590092", "EPR610827", "EPR617819", "EPR637054", "EPR649884", 
		"EPR776038", "EPR808798", "EPR833081", "EPR848660", "EPR892053", 
		"EPR910408", "EPR944574", "EPR024394", "EPR067643", "EPR085015", 
		"EPR109239", "EPR147621", "EPR152012", "EPR161836", "EPR189569", 
		"EPR222952", "EPR252773", "EPR269709", "EPR305820", "EPR321343", 
		"EPR375407", "EPR376132", "EPR414849", "EPR439346", "EPR439465", 
		"EPR441224", "EPR463296", "EPR475366", "EPR520858", "EPR529818", 
		"EPR535504", "EPR550639", "EPR568653", "EPR656841", "EPR733768", 
		"EPR737855", "EPR739383", "EPR739964", "EPR813110", "EPR821311", 
		"EPR843849", "EPR859838", "EPR869234", "EPR876626", "EPR932584", 
		"EPR006055", "EPR018999", "EPR048609", "EPR054525", "EPR055665", 
		"EPR087793", "EPR107163", "EPR121715", "EPR178912", "EPR188243",
		"EPR195073", "EPR200531", "EPR204507", "EPR229625", "EPR232859", 
		"EPR274814", "EPR275423", "EPR287449", "EPR318865", "EPR320033", 
		"EPR322524", "EPR356115", "EPR367652", "EPR367855", "EPR406645", 
		"EPR415517", "EPR438363", "EPR443196", "EPR451668", "EPR460913", 
		"EPR471772", "EPR484768", "EPR512820", "EPR545077", "EPR548213", 
		"EPR554926", "EPR562014", "EPR578337", "EPR590166", "EPR599137", 
		"EPR599897", "EPR626267", "EPR653593", "EPR670672", "EPR685442", 
		"EPR687508", "EPR687753", "EPR710404", "EPR719414", "EPR728793", 
		"EPR752649", "EPR757415", "EPR761353", "EPR774721", "EPR790594", 
		"EPR791050", "EPR807751", "EPR825019", "EPR829970", "EPR844841", 
		"EPR856363", "EPR864354", "EPR874182", "EPR884301", "EPR886748", 
		"EPR908949", "EPR909998", "EPR927666", "EPR945135", "EPR960265", 
		"EPR963250", "EPR966103", "EPR989096" ]

	print("LEN(to_process):", len(to_process))

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
	for item in to_process:
		clean_to_process.append({
			'file': 'hdfs://10.1.11.7:27000/test/int-plat/VCPLAT-2260/{}.{}.g.vcf.bgz'.format(item, chrom)
			'id': item,
			'pid': item
		})

	# Get version of sparse matrix
	version = path.basename(path.normpath(sparse_path))
	base = sparse_path.replace(version, '')
	self.log.debug('> Detected version of sparse matrix {}'.format(version))

	try:
		sm = hl.read_matrix_table(_name_with_chrom(sparse_path, chrom))
		self.log.info('> Sparse matrix {}/chrom-{} was loaded'.format(version, chrom))
		sm_loaded = True
	except:
		self.log.info('> Sparse matrix {}/chrom-{} could not be found and will be created'.format(version, chrom))
		sm_loaded = False

	import sys
	sys.exit()
	
	# Check for loaded experiments
	if sm_loaded:
		x = [ y.get('s') for y in self.data.col.collect() ]
		self.log.debug('> Loaded sparse matrix contains {} experiments'.format(len(x)))
		y = [ z for z in clean_to_process if z['id'] in x ]
		if len(y) != 0:
			self.log.error('> {} experiments are already loaded'.format(len(y)))
			clean_to_process = [ z for z in clean_to_process if z['id'] not in x ]

	# Create batches of samples to be loaded
	self.log.info('> Starting step 1 - creation of cumulative matrices of {} experiments, incrementing {} experiments at a time'.format(largeBatch, smallBatch))
	batches = _create__batches(clean_to_process, version, largeBatch, smallBatch)
		
	for idx1, batch in enumerate(batches):
		print('> Processing large batch {}/{} {}'.format(idx1, len(batches), batch[ 'version' ]))

		accum = None
		for idx2, pack in enumerate(batch[ 'content' ]):
			vsr = pack[ 'version' ]
			if idx2 == len(batch[ 'content' ]) - 1:
				vsr = batch[ 'version' ]
			small_batch_path = _name_with_chrom(path.join(base, vsr), chrom)
			print('     > Loading pack #{} of {} gVCF ({})'.format(idx2, len(batch[ 'content' ]), small_batch_path))
			for f in pack['content']:
				print(f)
			_load_gvcf(self.hl, pack[ 'content' ], small_batch_path, accum, chrom, config[ 'applications/combine/partitions_chromosome' ])
			accum = small_batch_path

	# Collect all the small sparse matrix and iteratively accumulate them
	revisions_to_collect = [ pack[ 'version' ] for pack in batches ]
	if sm_loaded:
		revisions_to_collect = [ (None, version) ] + revisions_to_collect

	print('>> SM TO COLLECT <<')
	#print(revisions_to_collect)

	self.log.info('> Starting step 2 - merging {} cumulative matrices'.format(len(revisions_to_collect)))
	for ii in range(1, len(revisions_to_collect)):
		print(ii, revisions_to_collect[ ii ])
		_combine_mt(self.hl, base, revisions_to_collect[ ii-1 ], revisions_to_collect[ ii ], utils.version_bump(revisions_to_collect[ ii ], 'version'), chrom)

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


def _load_gvcf(hl, experiments, version_path, previous_version_path, chrom, partitions):
	def transformFile(mt):
		return transform_gvcf(mt.annotate_rows(
			info = mt.info.annotate(MQ_DP = hl.null(hl.tint32), VarDP = hl.null(hl.tint32), QUALapprox = hl.null(hl.tint32))
		))
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
	comb.write(version_path, overwrite = True)


def _create__batches(experiments, version, largeSize = 500, smallSize = 100):
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