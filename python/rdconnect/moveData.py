import json
import requests

"""moveData

This module contains the functions used to move data from main cluster to HDFS.
"""

import os
from subprocess import call
from traceback import format_exc
from rdconnect.classLog import VoidLog

from rdconnect.utils import chrom_str_to_int, create_chrom_filename

def gvcf(config, log = VoidLog(), batch = 500):
	"""Function used to move (g)VCF files from a POSIX to HADOOP (HDFS) file 
	system.

	It makes a query to DM using the endpoint 'api_exp_status_list' to get the
	list of unprocessed experiments.

	Then, it wors as descrived in 'gvcf_from_file'. A summary follows:

	Using the values of 'moving_from' and 'moving_to' it creates the path to
	locate the gVCF to be moved. The real file to be moved from POSIX to HDSF,
	replacing as:

		[patient-id].chromosome.g.vcf.gz

	Where 'chromosome' will be replaced by the current chromosome being
	transferred.

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
	batch: number of files to be moved as a batch. Once this limit is reached 
		the function ends.
	"""
	chrom = config['process/chrom']
	chrom = chrom_str_to_int(chrom)
	source_path = config['process/moving_from'], 
	destination_path = config['process/moving_to'], 
	cmd_1 = config['process/moving_s1'], 
	cmd_2main_conf['process/moving_s2']

	url = config['applications/datamanagement/ip']
	if not url.startswith('http://') and not url.startswith('https://'):
		url = 'https://{0}'.format(url)

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
	data = "{\"" + flag + "\": \"" + value + "\"}"


	data = "{\"page\": 1, \"pageSize\": " + str(batch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"], \"sorted\":[{\"id\":\"RD_Connect_ID_Experiment\",\"desc\":false}], \"filtered\":[{\"id\":\"variantCalling\",\"value\":\"pass\"},{\"id\":\"rohs\",\"value\":\"pass\"},{\"id\":\"in_platform\",\"value\":\"waiting\"}]}"

	log.debug('> Querying DM using url "{1}"'.format(url))
	
	response = requests.post(url, data = data, headers = headers, verify = False)
	if response.status_code != 200:
		log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
	


def gvcf_from_file(config, log = VoidLog()):
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
	gvcf_list = config['process/moving_gvcf_list'], 
	source_path = config['process/moving_from'], 
	destination_path = config['process/moving_to'], 
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
				line = '{}.chromosome.g.vcf.gz'.frormat(line.strip().split('/')[ 1 ])
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

