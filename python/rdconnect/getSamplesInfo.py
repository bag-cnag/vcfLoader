
import json
import requests

from rdconnect.classLog import VoidLog

def experiments(config, log = VoidLog(), is_playground = False):
	"""Function used to query data-management and get a list of all experiments.

	The queering is done using the token provided by 
	'applications/datamanagement/token' to the end point provided by the 
	configuration slot:

		'applications/datamanagement/ip' ('applications/datamanagement/host')

	The base URL for the API consumption is indicated, using the argument 
	'is_playground', at:

		'applications/datamanagement/api_exp'
		'applications/datamanagement/api_exp_playground'

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job.
	log: logger, optional
		Used to track the queering process
	is_playground: bool, optional
		Boolean indicating it playground URLs shall be used.
	"""
	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = config['applications/datamanagement/api_exp_playground'].format(url)
	else:
		url = config['applications/datamanagement/api_exp'].format(url)
	headers = { 
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	log.debug('Querying experiment by group using url "{}"'.format(url))
	resp = requests.get(url, headers = headers, verify = False)
	data = json.loads(resp.content)
	return data


def experiment_by_group(config, log = VoidLog(), is_playground = False):
	"""Function used to query data-management and get a list of experiments
	that belongs to a specific group of users.

	The queering is done using the token provided by 
	'applications/datamanagement/token' to the end point provided by the 
	configuration slot:

		'applications/datamanagement/ip' ('applications/datamanagement/host')

	The base URL for the API consumption is indicated, using the argument 
	'is_playground', at:

		'applications/datamanagement/api_exp_group'
		'applications/datamanagement/api_exp_group_playground'

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job that must include the keys 
		'applications/datamanagement/ip','applications/datamanagement/token',
		'applications/datamanagement/host', 'applications/datamanagement/api_exp_group',
		and 'applications/datamanagement/api_exp_group_playground'.
	log: logger, optional
		Used to track the queering process
	is_playground: bool, optional
		Boolean indicating it playground URLs shall be used.
	"""
	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = config['applications/datamanagement/api_exp_group_playground'].format(
			url, config['applications/api_group'])
	else:
		url = config['applications/datamanagement/api_exp_group'].format(
			url, config['applications/api_group'])
	headers = { 
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	log.debug('Querying experiment by group using url "{}"'.format(url))
	resp = requests.get(url, headers = headers, verify = False)
	data = json.loads(resp.content)
	return data


def experiment_status(config, log = VoidLog(), is_playground = False):
	"""Function used to query data-management and get the status of a the 
	experiments belonging to a group of users.

	The queering is done using the token provided by 
	'applications/datamanagement/token' to the end point provided by the 
	configuration slot:

		'applications/datamanagement/ip' ('applications/datamanagement/host')

	The base URL for the API consumption is indicated, using the argument 
	'is_playground', at:

		'applications/datamanagement/api_exp_status'
		'applications/datamanagement/api_exp_status_playground'

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job that must include the keys 
		'applications/datamanagement/ip','applications/datamanagement/token',
		'applications/datamanagement/host', 'applications/datamanagement/api_exp_status',
		and 'applications/datamanagement/api_exp_status_playground'.
	log: logger, optional
		Used to track the queering process
	is_playground: bool, optional
		Boolean indicating it playground URLs shall be used.
	"""
	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = config['applications/datamanagement/api_exp_status_playground'].format(
			url, config['applications/api_group'])
	else:
		url = config['applications/datamanagement/api_exp_status'].format(
			url, config['applications/api_group'])
	headers = { 
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}	
	log.debug('Querying experiment\'s status using url "{}"'.format(url))
	resp = requests.get(url, headers = headers, verify = False)
	data = json.loads(resp.content)
	return data



def experiments_to_process(experiment_available, experiment_status, check_hdfs = False):
	"""Given the experiments seen by the user as well as their status, returns 
	the ones that are in HDFS and have to be processed.

	It looks for those experiments that had 'genomicsdb' set to 'waiting'.

	Parameters
	----------
	experiment_available: list, mandatory
		List created with the function 'experiment_by_group'
	experiment_status: list, mandatory
		List created with the function 'experiment_status'
	check_hdfs: bool, optional
		By default 'False'. If set to true if filters out those experiments
		where 'hdfs' is not set to 'pass'.
	"""
	experiment_status = [ x for x in experiment_status if x[ 'genomicsdb' ] == 'waiting' ]
	if check_hdfs:
		experiment_status = [ x for x in experiment_status if x[ 'hdfs' ] == 'pass' ]
	experiment_status_2 = [ x[ 'Experiment' ] for x in experiment_status ]
	experiment_available_2 = [ x[ 'RD_Connect_ID_Experiment' ] for x in experiment_available ]
	selected_experiments = [ x for x in experiment_available_2 if x in experiment_status_2 ]
	return experiment_available
	

def experiments_and_family(pids, config, sort_output = True):
	url = config['applications/phenostore/api_exp_mul'].format(config['applications/phenostore/ip'])
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)
	headers = { 'Content-Type': 'application/json', 
		'Authorization': config['application/kc_token'], 
		'Host': config['applications/phenostore/host'] 
	}
	data = []
	print('pids': pids)
	for ii in range(0,(len(pids)//1000)+1) :
		#body = { 'patients': [ { 'id': x[ 'Participant_ID' ] } for x in pids[(i*1000):((i+1)*1000)] ] }
		body = { 'patients': [ { 'id': x[ 1 ] } for x in pids[(ii*1000):((ii+1)*1000)] ] }
		print("--->", url)
		print("---with-->", body)
		resp = requests.post(url, headers = headers, json = body, verify = False)
		print(resp.text)
		x = resp.json()
		if x is dict:
			raise Exception("Communication with PhenosTore giving \"{}\" returned 'dict' instead of 'list' with content \"{}\"".format(str(body), str(x)))
		data.append(x)
	print(data)
	parsed = {}
	for elm in data:
		pid = list( elm.keys() )[ 0 ]
		if type( elm[ pid ] ) == str:
			fam = '---'
		else: 
			fam = elm[ pid ][ 'family' ] if 'family' in elm[ pid ].keys() else '---'
		parsed[ pid ] = fam
	rst = [ [ pak[ 0 ], pak[ 1 ], parsed[ pak[ 1 ] ] ] for pak in pids ]
	if sort_output:
		return sorted( rst, key = lambda x: x[ 2 ] )
	else:
		return rst




def experiments_with_dm_traking(pids, config, log, n = 200):
	packs = []
	for ii in range(0, len(pids), n):  
		packs.append(','.join(pids[ii:ii + n]))
	log.debug('> Data-management will be queried {} times, each time with {} experiments'.format(len(packs), n))

	url = config['applications/datamanagement/api_sm'].format(config['applications/datamanagement/ip'])
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	headers = { 'accept': 'application/json', 
		'Content-Type': 'application/json', 
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] }

	log.debug('> Created query URL for data-management: {}'.format(url))

	table = {}
	for ii, samlist in enumerate(packs):
		q_url = url + '?experiments=' + samlist
		response = requests.get(q_url, headers = headers, verify = False)
		if response.status_code != 200:
			log.error('> Data-management returned {} ("{}") when queried with #{} batch of experiments'.format(response.status_code, response.text, ii))
			return 
		else:
			data = json.loads(response.content)
			table.update(data)
	return table