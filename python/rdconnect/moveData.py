"""moveData

This module contains the functions used to move data from main cluster to HDFS.
"""

import os
from subprocess import call
from traceback import format_exc

from utils import chrom_str_to_int, create_chrom_filename

def move_gvcf(log, chrom, gvcf_list, source_path, destination_path, cmd_1, cmd_2):
	"""Function used to move (g)VCF files from a POSIX to HADOOP (HDFS) file 
	system.

	The content of the 'gvcf_list' corresponds to a fragment of the path
	used to locate the origin of the (g)VCF files. The content of this file
	must follow this schema:

		[group-name]/[patient-id]

	The real file to be moved from POSIX to HDSF will be looked at path 
	indicated by 'source_path', replacing the string 'filename' by:

		[patient-id].chromosome.g.vcf.gz

	Where 'chromosome' will be replaced by the current chromosome being
	transferred.

	This function relies on to commands ('cmd_1' and 'cmd_2') to indicate how
	to perform the multiple ssh commands. An example of them could be:

		1) "ssh -i /home/architect/.ssh/id_rsa trinity@10.2.2.2 ' cat ",
 		2) "ssh -i /home/architect/.ssh/id_rsa trinity@hdfs 'source /home/trinity/.bash_profile; hdfs dfs -put  -f - "			

	Parameters
	----------
	chrom: str, mandatory
		The chromosome for which (g)VCF files will be transferred
	gvcf_list: context, mandatory
		HAIL context
	source_path: str, mandatory
		List of groups and files to be moved
	source_path: str, mandatory
		String used to save the dense matrix
	destination_path: str, mandatory
		String used to save the dense matrix
	cmd_1: str, mandatory
		String used to save the dense matrix
	cmd_2: str, mandatory
		String used to save the dense matrix
	"""
	log.info('Entering step "move_gvcf"')
	log.debug('- Argument "chrom" filled with "{}"'.format(chrom))
	log.debug('- Argument "gvcf_list" filled with "{}"'.format(gvcf_list))
	log.debug('- Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('- Argument "destination_path" filled with "{}"'.format(destination_path))
	log.debug('- Argument "cmd_1" filled with "{}"'.format(cmd_1))
	log.debug('- Argument "cmd_2" filled with "{}"'.format(cmd_2))
	chrom = chrom_str_to_int(chrom)
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

