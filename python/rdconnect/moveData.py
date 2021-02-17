import os
import sys
import json
import requests
from os import path

import rdconnect.getSamplesInfo as get
from subprocess import call
from traceback import format_exc
from rdconnect.classLog import VoidLog
from rdconnect.utils import chrom_str_to_int, create_chrom_filename

"""moveData

This module contains the functions used to move data from main cluster to HDFS.
"""

def get_experiments_prepared(config, log = VoidLog(), batch = 500, is_playground = False):
	"""Function to create a list of experiments that have to be moved to HDFS or
	CEPH.

	This function calls the data manager to get a batch of files ready to be
	processed, starting to be moved to the distributed file system. The call
	relies on the following filtering criteria:

		variantCalling: pass
		hdfs: waiting
		genomicsdb: waiting
		es: waiting

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job that must include the keys 'process/chrom',
		'process/moving_from', 'process/moving_to',
		'process/moving_s1', and 'process/moving_s2'.
	log: logger, optional
		Used to track the moving process
	batch: numeric, optional
		Number of files to be moved as a batch. Once this limit is reached 
		the function ends (default: 500).
	
	Returns
	-------
	Nothing. The list of files is written to disk as "transfer_files.txt".
	"""
	if config is None:
		raise Exception('Started "get_experiments_prepared" and no "config" was provided.')

	chrm_str = config['process/chrom']
	if chrm_str is None:
		chrm_str = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'.split(',')
	else:
		chrm_str = chrm_str.split(',')
	source_path = config['process/moving_from']
	destination_hdfs = config['process/moving_to_hdfs']
	destination_ceph = config['process/moving_to_ceph']

	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = config['applications/datamanagement/api_exp_status_list_playground'].format(url)
	else:
		url = config['applications/datamanagement/api_exp_status_list'].format(url)

	log.info('Entering step "get_experiments_prepared"')
	log.debug('> Argument "chrom" filled with "{}"'.format(chrm_str))
	log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('> Argument "destination_hdfs" filled with "{}"'.format(destination_hdfs))
	log.debug('> Argument "destination_ceph" filled with "{}"'.format(destination_ceph))

	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	
	data = "{\"page\": 1, \"pageSize\": "  + str(batch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"],\"sorted\": [{\"id\": \"RD_Connect_ID_Experiment\",\"desc\": false}],\"filtered\": [{\"id\": \"variantCalling\",\"value\": \"pass\"},{\"id\": \"hdfs\",\"value\": \"waiting\"},{\"id\": \"genomicsdb\",\"value\": \"waiting\"},{\"id\": \"es\",\"value\": \"waiting\"}]}"
	log.debug('> Querying DM using URL "{0}"'.format(url))

	log.debug(headers)
	log.debug(data)
	log.debug(url)
	
	response = requests.post(url, data = data, headers = headers, verify = False)
	if response.status_code != 200:
		log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
		sys.exit(2)

	to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
	log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))
	
	all_group = get.experiment_by_group(config, log, is_playground)
	log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))
	
	to_process_group = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]
	with open('transfer_files.txt', 'w') as fw:
		for ii, xx in enumerate(to_process_group):
			for chrm in chrm_str:
				c1 = source_path.replace('[patient-id]', xx['RD_Connect_ID_Experiment']).replace('[owner]', xx['Owner']).replace('[chromosome]', chrm)
				c2 = destination_ceph.replace('[patient-id]', xx['RD_Connect_ID_Experiment']).replace('[owner]', xx['Owner']).replace('[chromosome]', chrm)
				c3 = destination_hdfs.replace('[patient-id]', xx['RD_Connect_ID_Experiment']).replace('[owner]', xx['Owner']).replace('[chromosome]', chrm)
				c4 = xx['RD_Connect_ID_Experiment']
				fw.write(c1 + '\t' + c2 + '\t' + c3 + '\t' + c4 + '\n')




def gvcf_to_hdfs(config, log = VoidLog(), batch = 500, include_tbi = True, is_playground = False):
	"""Function used to move (g)VCF files from a POSIX to HADOOP (HDFS) file 
	system.

	It makes a first query to DM using the endpoint 'api_exp_status_list' to 
	get the list of unprocessed experiments. The it used the function 
	'experiment_by_group' from the module 'getSamplesInfo' to get their owner.

	Then, and taking the sour path from 'moving_from' it will add the location
	of the file. The bracket-names in that variables will be substituted, as 
	well as the ones from 'moving_to', that indicates where the file will be
	places, as:

		* [owner]: name of the PI/Group/... owing the experiment
		* [patient-id]: identifier of the experiment/patient
		* [chromosome]: numeric value for the chromosome

	This function relies on to commands ('moving_s1' and 'moving_s2') to 
	indicate how to perform the multiple ssh commands. An example of them 
	could be:

		1) "ssh -i /home/architect/.ssh/id_rsa trinity@10.2.2.2 ' cat ",
 		2) "ssh -i /home/architect/.ssh/id_rsa trinity@hdfs 'source /home/trinity/.bash_profile; hdfs dfs -put  -f - "			

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job that must include the keys 'process/chrom',
		'process/moving_from', 'process/moving_to',
		'process/moving_s1', and 'process/moving_s2'.
	log: logger, optional
		Used to track the moving process
	batch: numeric, optional
		Number of files to be moved as a batch. Once this limit is reached 
		the function ends (default: 500).
	include_tbi, boolean, optional
		If set to True it also copies the tbi files associated to the gVCF 
		(default: True).
	Returns
	-------
	A list with the "RD_Connect_ID_Experiment" for moved experiments.
	"""
	chrom = chrom_str_to_int(str(config['process/chrom']))
	source_path = config['process/moving_from']
	destination_path = config['process/moving_to']
	cmd_1 = config['process/moving_s1']
	cmd_2 = config['process/moving_s2']

	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

	if is_playground:
		url = config['applications/datamanagement/api_exp_status_list_playground'].format(url)
	else:
		url = config['applications/datamanagement/api_exp_status_list'].format(url)

	log.info('Entering step "gvcf"')
	log.debug('> Argument "chrom" filled with "{}"'.format(chrom))
	log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	log.debug('> Argument "cmd_1" filled with "{}"'.format(cmd_1))
	log.debug('> Argument "cmd_2" filled with "{}"'.format(cmd_2))

	headers = { 
		'accept': 'application/json', 'Content-Type': 'application/json',
		'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
		'Host': config['applications/datamanagement/host'] 
	}
	data = "{\"page\": 1, \"pageSize\": " + str(batch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"], \"sorted\":[{\"id\":\"RD_Connect_ID_Experiment\",\"desc\":false}], \"filtered\":[{\"id\":\"variantCalling\",\"value\":\"pass\"},{\"id\":\"rohs\",\"value\":\"pass\"},{\"id\":\"in_platform\",\"value\":\"waiting\"}]}"
	log.debug('> Querying DM using url "{0}"'.format(url))

	log.debug(headers)
	log.debug(data)
	log.debug(url)
	
	response = requests.post(url, data = data, headers = headers, verify = False)
	if response.status_code != 200:
		log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
		sys.exit(2)

	to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
	log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))
	
	all_group = get.experiment_by_group(config, log, is_playground)
	log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))
	
	to_process_group = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]


	# chrom_str = chrom
	# if chrom_str == '23':
	# 	chrom_str = 'MT'
	# elif chrom_str == '24':
	# 	chrom_str = 'X'
	# elif chrom_str == '25':
	# 	chrom_str = 'Y'
	
	# for idx, line in enumerate(to_process_group):
	# 	log.debug('Processing samples #{} "{}"'.format(str(idx), line['RD_Connect_ID_Experiment']))
	# 	try:
	# 		ori = source_path.replace('[owner]', line['Owner'])\
	# 			.replace('[patient-id]', line['RD_Connect_ID_Experiment'])\
	# 			.replace('[chromosome]', str(chrom_str))
	# 		des = destination_path.replace('[owner]', line['Owner'])\
	# 			.replace('[patient-id]', line['RD_Connect_ID_Experiment'])\
	# 			.replace('[chromosome]', str(chrom))	

	# 		log.debug('>> Moving experiment {} from "{}" to "{}"'.format(line['RD_Connect_ID_Experiment'], ori, des))
			
	# 		command_1 = cmd_1 + ori + " ' | "
	# 		command_2 = cmd_2 + des + "'"
	# 		command = command_1 + command_2
	# 		os.system(command)

	# 		if include_tbi:
	# 			command_1 = cmd_1 + ori + ".tbi ' | "
	# 			command_2 = cmd_2 + des + ".tbi'"
	# 			command = command_1 + command_2
	# 			os.system(command)
	# 	except Exception as ex:
	# 		log.error('Unexpected error:\n{}'.format(str(ex)))
	# 		log.debug('Stack: {}'.format(str(format_exc())))
	# 		sys.exit(2)

	# return [ x['RD_Connect_ID_Experiment'] for x in to_process_group ]
	


def gvcf_to_hdfs_from_file(config, log = VoidLog()):
	"""Function used to move (g)VCF files from a POSIX to HADOOP (HDFS) file 
	system.

	The content of the 'moving_gvcf_list' corresponds to a fragment of the path
	used to locate the origin of the (g)VCF files. The content of this file
	must follow this schema:

		[group-name]/[patient-id]

	The real file to be moved from POSIX to HDSF will be looked at path 
	indicated by 'source_path', replacing the string 'filename' by:

		[patient-id].chromosome.g.vcf.gz

	Where 'chromosome' will be replaced by the current chromosome being
	transferred.

	This function relies on to commands ('moving_s1' and 'moving_s1') to 
	indicate how to perform the multiple ssh commands. An example of them 
	could be:

		1) "ssh -i /home/architect/.ssh/id_rsa trinity@10.2.2.2 ' cat ",
 		2) "ssh -i /home/architect/.ssh/id_rsa trinity@hdfs 'source /home/trinity/.bash_profile; hdfs dfs -put  -f - "			

	Parameters
	----------
	config: ConfigFile, mandatory
		Configuration for the job that must include the keys 'process/chrom',
		'process/moving_gvcf_list', 'process/moving_from', 'process/moving_to',
		'process/moving_s1', and 'process/moving_s2'.
	log: logger, optional
		Used to track the moving process
	"""
	chrom = config['process/chrom']
	chrom = chrom_str_to_int(chrom)
	gvcf_list = config['process/moving_gvcf_list']
	source_path = config['process/moving_from']
	destination_path = config['process/moving_to']
	cmd_1 = config['process/moving_s1'],
	cmd_2main_conf['process/moving_s2']

	log.info('Entering step "gvcf_from_file"')
	log.debug('> Argument "chrom" filled with "{}"'.format(chrom))
	log.debug('> Argument "gvcf_list" filled with "{}"'.format(gvcf_list))
	log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	log.debug('> Argument "cmd_1" filled with "{}"'.format(cmd_1))
	log.debug('> Argument "cmd_2" filled with "{}"'.format(cmd_2))
	
	with open(gvcf_list, 'r') as rd:
		for idx, line in enumerate(rd):
			log.debug('Processing line #{} with content "{}"'.format(str(idx), line.strip()))
			try:
				line = '{}.chromosome.g.vcf.gz'.format(line.strip().split('/')[ 1 ])
				line = create_chrom_filename(line, chrom)
				ori = source_path.replace('filename', line)
				des = os.path.join(destination_path, line).replace('gz', 'bgz')

				command_1 = cmd_1 + ori + " ' | "
				command_2 = cmd_2 + des + "'"
				command = command_1 + command_2
				os.system(command)
			except Exception as ex:
				log.error('Unexpected error:\n{}'.format(str(ex)))
				log.debug('Stack: {}'.format(str(format_exc())))
				sys.exit(2)

