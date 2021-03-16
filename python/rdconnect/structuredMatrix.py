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

	#experiments = "\"E708461\", \"E913959\", \"E581219\", \"E818684\", \"E252031\", \"E835057\", \"E502175\", \"E442076\", \"E477152\", \"E984397\", \"E685936\", \"E670179\", \"E776395\", \"E602052\", \"E677263\", \"E518789\", \"E361755\", \"E093681\", \"E374480\", \"E459670\", \"E819495\", \"E782614\", \"E001607\", \"E067269\", \"E564489\", \"E054978\", \"E518159\", \"E768938\", \"E277400\", \"E131378\", \"E200174\", \"E884282\", \"E716489\", \"E680922\", \"E466501\", \"E763199\", \"E659464\", \"E496283\", \"E398007\", \"E065995\", \"E652764\", \"E840290\", \"E202537\", \"E869976\", \"E241429\", \"E961966\", \"E683302\", \"E227425\", \"E138363\", \"E819415\", \"E823519\", \"E912622\", \"E228435\", \"E836596\", \"E736187\", \"E036579\", \"E182604\", \"E850617\", \"E364158\", \"E646578\", \"E335360\", \"E855503\", \"E941761\", \"E556858\", \"E226112\", \"E293199\", \"E484644\", \"E025768\", \"E394768\", \"E953392\", \"E529918\", \"E687071\", \"E569670\", \"E715220\", \"E250812\", \"E708296\", \"E342739\", \"E060385\", \"E632096\", \"E644156\", \"E896806\", \"E907380\", \"E763057\", \"E705220\", \"E587584\", \"E385883\", \"E156659\", \"E090202\", \"E002768\", \"E355526\", \"E908575\", \"E703715\", \"E113473\", \"E985246\", \"E656687\", \"E491299\", \"E950584\", \"E743766\", \"E430650\", \"E990345\", \"E451978\", \"E818146\", \"E015348\", \"E955106\", \"E637633\", \"E862785\", \"E503981\", \"E831395\", \"E018706\", \"E842462\", \"E196012\", \"E546623\", \"E111335\", \"E385255\", \"E201269\", \"E124628\", \"E918922\", \"E378336\", \"E387165\", \"E019027\", \"E142452\", \"E516063\", \"E013668\", \"E536145\", \"E802988\", \"E060486\", \"E126548\", \"E825508\", \"E316441\", \"E524007\", \"E295108\", \"E759739\", \"E442820\", \"E109017\", \"E351371\", \"E546195\", \"E204857\", \"E974965\", \"E372193\", \"E575431\", \"E146992\", \"E076648\", \"E863352\", \"E610329\", \"E582574\", \"E496220\", \"E050920\", \"E341411\", \"E427508\", \"E028959\", \"E504859\", \"E318804\", \"E254625\", \"E785419\", \"E465415\", \"E233531\", \"E935151\", \"E927008\", \"E748847\", \"E530654\", \"E485203\", \"E268393\", \"E630047\", \"E082430\", \"E873977\", \"E711971\", \"E419757\", \"E395178\", \"E750770\", \"E406997\", \"E568337\", \"E547969\", \"E500206\", \"E145586\", \"E972990\", \"E585081\", \"E053668\", \"E129927\", \"E831461\", \"E017496\", \"E697292\", \"E607221\", \"E172432\", \"E055675\", \"E242338\", \"E295663\", \"E088110\", \"E420087\", \"E084538\", \"E919182\", \"E997098\", \"E182182\", \"E359934\", \"E588408\", \"E674448\", \"E679279\", \"E528144\", \"E373322\", \"E213047\", \"E900312\", \"E396592\", \"E694018\", \"E633365\", \"E630260\", \"E580632\", \"E641781\", \"E618374\", \"E314739\", \"E091846\", \"E606135\", \"E998401\", \"E307676\", \"E236809\", \"E328808\", \"E823118\", \"E401011\", \"E568909\", \"E597228\", \"E363045\", \"E634759\", \"E529001\", \"E180311\", \"E544587\", \"E369022\", \"E114126\", \"E722731\", \"E701278\", \"E223966\", \"E352382\", \"E077345\", \"E048964\", \"E067968\", \"E461348\", \"E733039\", \"E236206\", \"E653452\", \"E749305\", \"E961641\", \"E295343\", \"E264595\", \"E847007\", \"E231151\", \"E904864\", \"E341757\", \"E573252\", \"E301208\", \"E392833\", \"E551753\", \"E721249\", \"E955160\", \"E140117\", \"E021559\", \"E990493\", \"E629395\", \"E691845\", \"E677261\", \"E643847\", \"E540717\", \"E226591\", \"E135345\", \"E119370\", \"E230102\", \"E269050\", \"E006533\", \"E122002\", \"E205713\", \"E932287\", \"E342041\", \"E608630\", \"E328747\", \"E215088\", \"E024915\", \"E015811\", \"E458394\", \"E620926\", \"E666576\", \"E655845\", \"E719768\", \"E617775\", \"E198840\", \"E480759\", \"E375913\", \"E649888\", \"E727484\", \"E415657\", \"E764967\", \"E634933\", \"E833949\", \"E419012\", \"E785199\", \"E176886\", \"E374852\", \"E333740\", \"E222045\", \"E068972\", \"E828350\", \"E767398\", \"E394250\", \"E834997\", \"E834305\", \"E684691\", \"E162665\", \"E221158\", \"E178234\", \"E181736\", \"E026657\", \"E619208\", \"E515074\", \"E728960\", \"E323137\", \"E709240\", \"E593666\", \"E679757\", \"E344726\", \"E418152\", \"E790155\", \"E837346\", \"E857745\", \"E931031\", \"E815038\", \"E278971\", \"E530263\", \"E231174\", \"E346814\", \"E529993\", \"E962782\", \"E131811\", \"E057396\", \"E151042\", \"E263772\", \"E944255\", \"E593976\", \"E398014\", \"E199100\", \"E793671\", \"E900072\", \"E110124\", \"E972159\", \"E307532\", \"E868820\", \"E683753\", \"E035883\", \"E210112\", \"E536199\", \"E280915\", \"E259459\", \"E994014\", \"E660702\", \"E259439\", \"E832336\", \"E957794\", \"E576370\", \"E641168\", \"E852443\", \"E465935\", \"E845313\", \"E474740\", \"E399811\", \"E566703\", \"E770981\", \"E569714\", \"E061486\", \"E740617\", \"E344463\", \"E108679\", \"E356271\", \"E024963\", \"E650101\", \"E498215\", \"E096909\", \"E595072\", \"E069915\", \"E411236\", \"E747587\", \"E698340\", \"E234187\", \"E369353\", \"E394334\", \"E072637\", \"E863248\", \"E867196\", \"E781848\", \"E355248\", \"E672316\", \"E688261\", \"E146101\", \"E290233\", \"E490357\", \"E869829\", \"E245218\", \"E033124\", \"E670673\", \"E421038\", \"E435580\", \"E718485\", \"E501880\", \"E873255\", \"E809044\", \"E922440\", \"E909559\", \"E859035\", \"E755872\", \"E929787\", \"E601879\", \"E360678\", \"E368349\", \"E882904\", \"E907423\", \"E791014\", \"E072078\", \"E560737\", \"E225708\", \"E226686\", \"E994358\", \"E858200\", \"E563541\", \"E650957\", \"E172199\", \"E339464\", \"E575404\", \"E524074\", \"E712166\", \"E001307\", \"E188349\", \"E324155\", \"E320724\", \"E492336\", \"E044266\", \"E755784\", \"E053247\", \"E971029\", \"E411654\", \"E995692\", \"E727235\", \"E019777\", \"E990989\", \"E315039\", \"E463716\", \"E851513\", \"E039407\", \"E857778\", \"E859687\", \"E929105\", \"E960186\", \"E521504\", \"E368583\", \"E230026\", \"E170192\", \"E422490\", \"E123547\", \"E280445\", \"E049839\", \"E762160\", \"E683349\", \"E517619\", \"E419164\", \"E371797\", \"E378284\", \"E665005\", \"E518982\", \"E595437\", \"E836940\", \"E668836\", \"E338155\", \"E981846\", \"E670099\", \"E429812\", \"E302839\", \"E056115\", \"E610270\", \"E762641\", \"E824517\", \"E580618\", \"E794260\", \"E714719\", \"E459775\", \"E773605\", \"E147597\", \"E806203\", \"E668139\", \"E609481\", \"E247495\", \"E650968\", \"E604209\", \"E468656\", \"E034911\", \"E451189\", \"E579171\", \"E451711\", \"E029478\", \"E555135\", \"E490724\", \"E094157\", \"E612192\", \"E432355\", \"E606250\", \"E038410\", \"E282097\", \"E407955\", \"E533486\", \"E125450\", \"E858544\", \"E257079\", \"E154877\", \"E129682\", \"E458789\", \"E562030\", \"E138062\", \"E680215\", \"E722115\", \"E725538\", \"E203326\", \"E565759\", \"E969297\", \"E896100\", \"E658264\", \"E955002\", \"E233549\", \"E346722\", \"E173836\", \"E739287\", \"E268306\", \"E011770\", \"E166876\", \"E544922\", \"E900168\", \"E513847\", \"E106405\", \"E989225\", \"E648402\", \"E801718\", \"E112105\", \"E818091\", \"E155034\", \"E542178\", \"E348313\", \"E788510\", \"E273108\", \"E460876\", \"E157685\", \"E489834\", \"E233498\", \"E768818\", \"E935824\", \"E318856\", \"E765194\", \"E526392\", \"E794496\", \"E067227\", \"E341017\", \"E870081\", \"E239570\", \"E295013\", \"E250691\", \"E730903\", \"E424650\", \"E318312\", \"E639705\", \"E308470\", \"E856761\", \"E679475\", \"E670053\", \"E502433\", \"E169876\", \"E768467\", \"E096461\", \"E043987\", \"E830423\", \"E823080\", \"E948117\", \"E541794\", \"E290622\", \"E096528\", \"E190872\", \"E238943\", \"E183446\", \"E291718\", \"E296736\", \"E763666\", \"E809455\", \"E914376\", \"E011878\", \"E147906\", \"E073712\", \"E789294\", \"E859662\", \"E402757\", \"E939348\", \"E879313\", \"E762614\", \"E882743\", \"E361011\", \"E531423\", \"E369755\", \"E819011\", \"E688157\", \"E988638\", \"E271642\", \"E268937\", \"E885449\", \"E479947\", \"E018018\", \"E066925\", \"E180106\", \"E423909\", \"E133047\", \"E404917\", \"E453577\", \"E192035\", \"E276521\", \"E942263\", \"E343241\", \"E990404\", \"E554199\", \"E550295\", \"E033234\", \"E965145\", \"E004105\", \"E555848\", \"E942628\", \"E148111\", \"E358652\", \"E215194\", \"E566244\", \"E391479\", \"E386511\", \"E242167\", \"E798903\", \"E068769\", \"E019334\", \"E253322\", \"E784578\", \"E896291\", \"E177237\", \"E079677\", \"E524791\", \"E519421\", \"E639111\", \"E505044\", \"E778055\", \"E277574\", \"E551467\", \"E856322\", \"E361273\", \"E337032\", \"E602621\", \"E991668\", \"E091311\", \"E862472\", \"E238344\", \"E282498\", \"E442339\", \"E030378\", \"E561797\", \"E788354\", \"E507475\", \"E465545\", \"E652533\", \"E381147\", \"E585386\", \"E429491\", \"E035079\", \"E147688\", \"E762732\", \"E920434\", \"E080768\", \"E378799\", \"E263515\", \"E887912\", \"E066192\", \"E013555\", \"E369602\", \"E373936\", \"E641495\", \"E329801\", \"E510474\", \"E164263\", \"E142727\", \"E376302\", \"E190858\", \"E806188\", \"E560763\", \"E319814\", \"E121731\", \"E374387\", \"E434236\", \"E377552\", \"E821505\"".replace('"', '').split(",")
	#clean_to_process = [ { 'file': 'hdfs://10.1.11.7:27000/test/rdconnect/2717/{}.{}.g.vcf.bgz'.format(x.strip(), chrom_str), 'id': x.strip() } for x in experiments ]

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

	interval = utils.get_chrom_intervals(chrom, partitions, hl)
	vcfs = [ transformFile(mt) for mt in importFiles([ x[ 'file' ] for x in experiments ]) ]

	if previous_version_path == None:
		comb = combine_gvcfs(vcfs)
	else:
		previous = hl.read_matrix_table(previous_version_path)
		previous = previous.key_rows_by('locus')
		comb = combine_gvcfs([ previous ] + vcfs)

	#comb = comb.key_rows_by('locus', 'alleles')
	comb.write(version_path, overwrite = True)
	return comb


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
	dm_to_create = [ '0' ]
	print("to_add:", to_add)
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

