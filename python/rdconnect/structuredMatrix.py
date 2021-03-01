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
	experiments = ["E266728", 
		"E695592", 
		"E200125", 
		"E126992", 
		"E994789", 
		"E450955", 
		"E168211", 
		"E361148", 
		"E556777", 
		"E721161", 
		"E681340", 
		"E624770", 
		"E569305", 
		"E409646", 
		"E458620", 
		"E038082", 
		"E418033", 
		"E268295", 
		"E798212", 
		"E102231", 
		"E465324", 
		"E639597", 
		"E391859", 
		"E539071", 
		"E199064", 
		"E224456", 
		"E975451", 
		"E555485", 
		"E421857", 
		"E760631", 
		"E475295", 
		"E833307", 
		"E513714", 
		"E891378", 
		"E656412", 
		"E013667", 
		"E167937", 
		"E564491", 
		"E838803", 
		"E852852", 
		"E499287", 
		"E023124", 
		"E468551", 
		"E940848", 
		"E978419", 
		"E216835", 
		"E976219", 
		"E095595", 
		"E013810", 
		"E547462", 
		"E665526", 
		"E649555", 
		"E697785", 
		"E391396", 
		"E701314", 
		"E407103", 
		"E732369", 
		"E110248", 
		"E555470", 
		"E535275", 
		"E828989", 
		"E990887", 
		"E132073", 
		"E398838", 
		"E807156", 
		"E622571", 
		"E938693", 
		"E986975", 
		"E524798", 
		"E875496", 
		"E901505", 
		"E458866", 
		"E717835", 
		"E213123", 
		"E752474", 
		"E891911", 
		"E536984", 
		"E595785", 
		"E795741", 
		"E793493", 
		"E166593", 
		"E780021", 
		"E184978", 
		"E655971", 
		"E353725", 
		"E714666", 
		"E046090", 
		"E816335", 
		"E640640", 
		"E428318", 
		"E693681", 
		"E065787", 
		"E916101", 
		"E011892", 
		"E546599", 
		"E638904", 
		"E492844", 
		"E276404", 
		"E756790", 
		"E643318", 
		"E951644", 
		"E065444", 
		"E023742", 
		"E539389", 
		"E713810", 
		"E856191", 
		"E133853", 
		"E256145", 
		"E853071", 
		"E238807", 
		"E483708", 
		"E447484", 
		"E445814", 
		"E928419", 
		"E549240", 
		"E686717", 
		"E303680", 
		"E103285", 
		"E453103", 
		"E648846", 
		"E707461", 
		"E718607", 
		"E066469", 
		"E546671", 
		"E913046", 
		"E593387", 
		"E743504", 
		"E639042", 
		"E702891", 
		"E766258", 
		"E204573", 
		"E404450", 
		"E377072", 
		"E211313", 
		"E638901", 
		"E567058", 
		"E322164", 
		"E023116", 
		"E615073", 
		"E877195", 
		"E667934", 
		"E382378", 
		"E724485", 
		"E746298", 
		"E235233", 
		"E150791", 
		"E580380", 
		"E815826", 
		"E974162", 
		"E141859", 
		"E642566", 
		"E977199", 
		"E333091", 
		"E683023", 
		"E540661", 
		"E307836", 
		"E814673", 
		"E295021", 
		"E022269", 
		"E336273", 
		"E322150", 
		"E217306", 
		"E946844", 
		"E075972", 
		"E080164", 
		"E301086", 
		"E258342", 
		"E283773", 
		"E121986", 
		"E927494", 
		"E781049", 
		"E677231", 
		"E189842", 
		"E423948", 
		"E410174", 
		"E065452", 
		"E601526", 
		"E169735", 
		"E396618", 
		"E722419", 
		"E217569", 
		"E354978", 
		"E865859", 
		"E623449", 
		"E259608", 
		"E562813", 
		"E681028", 
		"E439851", 
		"E342635", 
		"E441144", 
		"E986336", 
		"E873608", 
		"E165572", 
		"E023638", 
		"E380036", 
		"E485372", 
		"E814740", 
		"E844816", 
		"E152371", 
		"E954695", 
		"E248166", 
		"E849426", 
		"E284043", 
		"E560215", 
		"E208118", 
		"E748563", 
		"E777305", 
		"E077316", 
		"E585338", 
		"E302358", 
		"E697145", 
		"E697767", 
		"E068244", 
		"E548363", 
		"E135984", 
		"E120426", 
		"E044567", 
		"E937473", 
		"E824319", 
		"E391888", 
		"E119169", 
		"E132020", 
		"E167930", 
		"E443968", 
		"E162533", 
		"E094259", 
		"E687886", 
		"E503061", 
		"E626322", 
		"E646146", 
		"E482371", 
		"E232386", 
		"E244568", 
		"E747743", 
		"E206118", 
		"E375130", 
		"E291671", 
		"E162434", 
		"E487331", 
		"E856671", 
		"E242745", 
		"E562627", 
		"E684555", 
		"E389612", 
		"E527643", 
		"E621578", 
		"E895023", 
		"E654350", 
		"E775603", 
		"E456355", 
		"E243615", 
		"E707003", 
		"E866183", 
		"E260502", 
		"E644779", 
		"E470502", 
		"E062265", 
		"E854394", 
		"E689254", 
		"E966709", 
		"E504027", 
		"E105418", 
		"E387510", 
		"E050099", 
		"E870808", 
		"E999730", 
		"E095055", 
		"E853217", 
		"E754827", 
		"E313065", 
		"E184644", 
		"E607955", 
		"E059864", 
		"E439866", 
		"E622209", 
		"E843237", 
		"E082356", 
		"E616325", 
		"E222035", 
		"E962634", 
		"E843655", 
		"E838657", 
		"E569362", 
		"E858943", 
		"E744369", 
		"E563976", 
		"E818177", 
		"E178988", 
		"E048102", 
		"E189356", 
		"E914217", 
		"E827321", 
		"E182944", 
		"E885371", 
		"E637864", 
		"E779446", 
		"E314629", 
		"E255279", 
		"E781931", 
		"E552938", 
		"E554301", 
		"E641234", 
		"E941429", 
		"E040007", 
		"E398435", 
		"E061608", 
		"E072965", 
		"E481417", 
		"E020374", 
		"E195585", 
		"E901884", 
		"E368021", 
		"E263189", 
		"E132117", 
		"E939034", 
		"E463437", 
		"E672247", 
		"E034079", 
		"E242645", 
		"E429315", 
		"E304340", 
		"E473086", 
		"E392377", 
		"E130814", 
		"E045518", 
		"E676905", 
		"E569770", 
		"E023482", 
		"E770282", 
		"E254578", 
		"E194636", 
		"E194003", 
		"E177144", 
		"E531898", 
		"E061560", 
		"E749941", 
		"E316487", 
		"E798422", 
		"E571257", 
		"E082804", 
		"E287761", 
		"E417197", 
		"E425194", 
		"E396195", 
		"E533847", 
		"E668211", 
		"E161638", 
		"E319118", 
		"E235794", 
		"E545998", 
		"E540519", 
		"E107975", 
		"E264808", 
		"E795585", 
		"E247456", 
		"E356754", 
		"E920716", 
		"E108836", 
		"E018547", 
		"E082370", 
		"E536842", 
		"E151669", 
		"E637184", 
		"E759708", 
		"E524271", 
		"E031095", 
		"E811279", 
		"E003447", 
		"E523670", 
		"E356577", 
		"E110634", 
		"E598292", 
		"E940135", 
		"E017735", 
		"E881245", 
		"E169188", 
		"E073507", 
		"E036657", 
		"E413016", 
		"E644262", 
		"E851222", 
		"E977246", 
		"E253723", 
		"E263222", 
		"E550012", 
		"E216977", 
		"E351450", 
		"E200452", 
		"E022991", 
		"E792716", 
		"E243669", 
		"E767712", 
		"E392768", 
		"E463408", 
		"E600861", 
		"E353137", 
		"E589284", 
		"E111081", 
		"E293514", 
		"E849549", 
		"E269678", 
		"E128884", 
		"E390070", 
		"E244628", 
		"E330932", 
		"E699677", 
		"E731490", 
		"E742274", 
		"E093068", 
		"E927475", 
		"E449359", 
		"E347325", 
		"E043673", 
		"E171070", 
		"E135072", 
		"E784471", 
		"E399844", 
		"E531794", 
		"E788429", 
		"E099660", 
		"E414102", 
		"E175420", 
		"E315915", 
		"E189666", 
		"E625652", 
		"E770533", 
		"E011914", 
		"E055718", 
		"E011821", 
		"E229656", 
		"E664689", 
		"E084338", 
		"E302529", 
		"E214365", 
		"E974393", 
		"E815069", 
		"E374724", 
		"E984972", 
		"E041780", 
		"E151410", 
		"E285516", 
		"E882258", 
		"E788819", 
		"E949500", 
		"E664232", 
		"E044978", 
		"E891655", 
		"E291831", 
		"E890091", 
		"E997743", 
		"E269306", 
		"E336458", 
		"E923722", 
		"E122162", 
		"E300966", 
		"E187639", 
		"E978742", 
		"E769556", 
		"E614182", 
		"E303887", 
		"E687156", 
		"E644271", 
		"E321169", 
		"E122205", 
		"E887674", 
		"E248679", 
		"E407725", 
		"E285034", 
		"E110296", 
		"E731325", 
		"E791049", 
		"E735017", 
		"E794645", 
		"E705064", 
		"E366522", 
		"E706285", 
		"E116814", 
		"E022624", 
		"E997066", 
		"E988891", 
		"E406586", 
		"E260302", 
		"E847680", 
		"E092597", 
		"E192028", 
		"E812074", 
		"E262081", 
		"E185575", 
		"E093809", 
		"E824863", 
		"E333676", 
		"E150300", 
		"E278087", 
		"E118176", 
		"E407604", 
		"E654208", 
		"E355807", 
		"E166619", 
		"E560878", 
		"E572505", 
		"E738600", 
		"E292484", 
		"E627226", 
		"E135279", 
		"E756110", 
		"E192088", 
		"E646472", 
		"E642143", 
		"E979451", 
		"E103641", 
		"E059377", 
		"E397013", 
		"E551175", 
		"E742729", 
		"E732471", 
		"E222347", 
		"E529491", 
		"E581929", 
		"E008457", 
		"E531124", 
		"E582283", 
		"E371394", 
		"E343918", 
		"E885226", 
		"E510971", 
		"E567343", 
		"E156237", 
		"E911338", 
		"E365102", 
		"E722809", 
		"E254166", 
		"E907751", 
		"E303328", 
		"E538437", 
		"E895554", 
		"E042557", 
		"E670973", 
		"E214606", 
		"E144483", 
		"E460586", 
		"E814054", 
		"E585439", 
		"E480733", 
		"E701756", 
		"E270715", 
		"E166634", 
		"E595981", 
		"E636750", 
		"E870304", 
		"E478568", 
		"E250925", 
		"E582714", 
		"E143221", 
		"E021333", 
		"E435170", 
		"E301318", 
		"E995796", 
		"E082594", 
		"E746705", 
		"E888099", 
		"E183749", 
		"E440909", 
		"E865570", 
		"E724551", 
		"E699923", 
		"E958394", 
		"E149675", 
		"E370952", 
		"E752078", 
		"E207024", 
		"E207024", 
		"E334753", 
		"E183433", 
		"E526883", 
		"E336575", 
		"E204203", 
		"E680920", 
		"E941600", 
		"E412038", 
		"E416431", 
		"E295185", 
		"E560176", 
		"E492532", 
		"E865773", 
		"E761279", 
		"E520873", 
		"E606154", 
		"E174651", 
		"E297642", 
		"E380476", 
		"E579799", 
		"E053220", 
		"E066580", 
		"E490493", 
		"E521488", 
		"E513607", 
		"E099037", 
		"E845333", 
		"E966984", 
		"E494609", 
		"E522316", 
		"E257268", 
		"E817104", 
		"E374199", 
		"E594658", 
		"E255090", 
		"E229697", 
		"E688153", 
		"E394224", 
		"E692030", 
		"E001090", 
		"E959610", 
		"E720659", 
		"E645599", 
		"E466428", 
		"E020513", 
		"E037866", 
		"E851997", 
		"E916529", 
		"E123123", 
		"E085416", 
		"E906171", 
		"E488675", 
		"E106820", 
		"E362648", 
		"E528989", 
		"E159781", 
		"E672047", 
		"E099941", 
		"E674374", 
		"E171459", 
		"E903390", 
		"E547627", 
		"E066274", 
		"E341465", 
		"E207938", 
		"E769498", 
		"E434138", 
		"E260914", 
		"E244392", 
		"E677072", 
		"E524824", 
		"E658454", 
		"E483974", 
		"E194925", 
		"E519013", 
		"E354439", 
		"E969610", 
		"E418776", 
		"E085257", 
		"E779523", 
		"E684358", 
		"E063665", 
		"E651208", 
		"E783057", 
		"E340565", 
		"E418096", 
		"E530091", 
		"E195170", 
		"E655211", 
		"E369244", 
		"E039786", 
		"E545599", 
		"E735890", 
		"E051214", 
		"E967317", 
		"E379763", 
		"E069306", 
		"E920776", 
		"E932090", 
		"E539329", 
		"E421160", 
		"E832060", 
		"E181803", 
		"E933003", 
		"E011121", 
		"E565312", 
		"E565392", 
		"E127955", 
		"E105565", 
		"E204603", 
		"E059160", 
		"E750555", 
		"E361309", 
		"E774529", 
		"E529099",
		"E455362", 
		"E663377", 
		"E361543", 
		"E423144", 
		"E249232", 
		"E110997", 
		"E899269", 
		"E016779", 
		"E595829", 
		"E399109", 
		"E940372", 
		"E953583", 
		"E509318", 
		"E423062", 
		"E954590", 
		"E802456", 
		"E972814", 
		"E266040", 
		"E449612", 
		"E639225", 
		"E202876", 
		"E155245", 
		"E502007", 
		"E582173", 
		"E686902", 
		"E080077", 
		"E854492", 
		"E029927", 
		"E477070", 
		"E742342", 
		"E702627", 
		"E542723", 
		"E840005", 
		"E853352", 
		"E248193", 
		"E914406", 
		"E153099", 
		"E482343", 
		"E178894", 
		"E296461", 
		"E321821", 
		"E220410", 
		"E434376", 
		"E725740", 
		"E140651", 
		"E079846", 
		"E494367", 
		"E463923", 
		"E671051", 
		"E674402", 
		"E786724", 
		"E480153", 
		"E548706", 
		"E828636", 
		"E719124", 
		"E468902", 
		"E636409", 
		"E585435", 
		"E887920", 
		"E106639", 
		"E339100", 
		"E638013", 
		"E792084" ]
	clean_to_process = []
	for item in experiments:
		if filesystem == 'ceph':
			clean_to_process.append({
				'file': 's3a://cnag/' + item[1],
				'id': item[3]
			})
		else:
			clean_to_process.append({
				#'file': item[2],
				#'id': item[3]
				'file': '/test/rdconnect/2545/{}.{}.g.vcf.bgz'.format(item, chrom),
				'id': item
			})

	# Get version of sparse matrix
	version = path.basename(path.normpath(sparse_path))
	base = sparse_path.replace(version, '')
	self.log.debug('> Detected version of sparse matrix {}'.format(version))

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
	
	for ii in range(1, len(revisions_to_collect)):
		last = _combine_mt(self.hl, base, revisions_to_collect[ ii-1 ], revisions_to_collect[ ii ], utils.version_bump(revisions_to_collect[ ii ][ 0 ], 'version'), chrom)

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
		comb = combine_gvcfs([ previous ] + vcfs)

	comb = comb.key_rows_by('locus', 'alleles')
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
	self.log.debug('Total of {0} experiments where read from file'.format(len(experiments)))
	print(experiments)

	idx = 0
	try:
		#for idx, batch in enumerate( mapping ):
		#	self.log.debug( "Flatting and filtering dense matrix {0} (sz: {1}) --> {2} - {3}".format( idx, len( batch ), batch[0], batch[len(batch) - 1] ) )
		#	sam = hl.literal([ x[ 0 ] for x in batch ], 'array<str>')
		sam = hl.literal(experiments_in_matrix, 'array<str>')
		small_matrix = sparse_matrix.filter_cols(sam.contains(sparse_matrix[ 's' ]))
		small_matrix = hl.experimental.sparse_split_multi(small_matrix, filter_changed_loci = True)
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



def dense_matrix_grouping(self = None, config = None, hl = None, log = VoidLog(), experiments = []):
	# self, isConfig, isHl = utils.check_class_and_config(None, config, hl, log, class_to = SparseMatrix)
	# self.log.info('Entering step "dense_matrix_grouping"')

	# if not isConfig:
	# 	self.log.error('No configuration was provided')
	# 	raise NoConfigurationException('No configuration was provided')

	# if not isHl:
	# 	self.log.error('No pointer to HAIL module was provided')
	# 	raise NoHailContextException('No pointer to HAIL module was provided')

	# self.log.debug('OVERWRITING chrom to chrom-21')

	# chrom = 21
	# sparse_path = path.join(self.config['applications/combine/sparse_matrix_path'], 'chrom-{}'.format(chrom))
	# sparse_matrix = hl.read_matrix_table(sparse_path)

	# experiments_in_matrix = [ x.get( 's' ) for x in sparse_matrix.col.collect() ]
	# self.log.debug('Obtained a total of {} experiments from sparse matrix (chrom 21)'.format(len(experiments_in_matrix)))

	# all_group = get.experiment_by_group(config, self.log, False)
	# self.log.debug('Obtained a total of {} experiments for the group'.format(len(all_group)))

	# full_ids_in_matrix = [ x for x in all_group if x[ 'RD_Connect_ID_Experiment' ] in experiments_in_matrix ]
	# print('full_ids_in_matrix', len( full_ids_in_matrix ))
	# print('\t', full_ids_in_matrix[ : 10 ])

	# experiments_and_families = get.experiments_and_family(full_ids_in_matrix, self.config)
	# print('experiments_and_families', len( experiments_and_families ))
	# print('\t', experiments_and_families[ : 10 ])

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

	full_samples = [ y.get('s') for y in self.data.col.collect() ]
	self.log.debug('> Number of samples in sparse matrix: {}'.format(len(full_samples)))
	self.log.debug('> First and last sample: {} // {}'.format(full_samples[0], full_samples[len(full_samples) - 1]))

	# packs = []
	# n = 200
	# for ii in range(0, len(full_samples), n):  
	# 	packs.append(','.join(full_samples[ii:ii + n]))
	# self.log.debug('> Data-management will be queried {} times, each time with {} experiments'.format(len(packs), n))

	# url = 'https://' + self.config['applications/datamanagement/api_sm'].format(self.config['applications/datamanagement/ip'])
	# headers = { #'accept': 'application/json', 
	# 	#'Content-Type': 'application/json', 
	# 	'Authorization': 'Token {0}'.format(self.config['applications/datamanagement/token']),
	# 	'Host': self.config['applications/datamanagement/host'] }

	# self.log.debug('> Created query URL for data-management: {}'.format(url))

	# table = {}
	# for ii, samlist in enumerate(packs):
	# 	q_url = url + '?experiments=' + samlist
	# 	response = requests.get(q_url, headers = headers, verify = False)
	# 	if response.status_code != 200:
	# 		self.log.error('> Data-management returned {} ("{}") when queried with #{} batch of experiments'.format(response.status_code, response.text, ii))
	# 		return 
	# 	else:
	# 		data = json.loads(response.content)
	# 		table.update(data)
	

	exp_dm = get.experiments_with_dm_traking(full_samples, self.config, self.log)
	print(exp_dm)
	exp_fam = get.experiments_and_family(full_samples, self.config)

	print(table)
	print(exp_fam)

	return self