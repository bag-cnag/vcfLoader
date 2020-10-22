"""loadGenetics

This module contains the functions used to load genetic data including:
	* Dense matrices of variants created with the functions in module 'moveData'
	* Germline VCF and gVCF files moved to HDFS with the functions in module 'moveData'
	* Somatic VCF and gVCF files moved to HDFS with the functions in module 'moveData'
	* CNV data stored as '.csv' files from HDFS, placed there using 'moveData' module functions
"""

import os
from traceback import format_exc

import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData
from rdconnect.classLog import VoidLog
from rdconnect.utils import _check_class_and_config

# #def dense_matrix(log, hl, source_path, destination_path = None):
# def dense_matrix(config, log, hl):
# 	"""This function is used to load a dense matrix into HAIL structure

# 	This function loads a dense matrix into HAIL structure in the HDFS. It 
# 	creates annotates the variants with general annotations:

# 		* samples, ad, dp, gtInt, gt, and gq freqInt, pos, alt, ref)

# 	Ant is also applied some filtering on the variants:

# 		* row filtering: both genotypes are reference, dp > 10, and gq > 20
# 		* remove duplicated according to locus and allele

# 	Parameters
# 	----------
# 	config: ConfigFile, mandatory
# 		Object of class 'ConfigFile' with the configuration of the pipeline. It
# 		will be used to extract 'source_path' and 'destination_path',
# 	log: logger, mandatory
# 		A logger to have track of the steps used in the loading process
# 	hl: context, mandatory
# 		HAIL context

# 	Returns
# 	-------
# 	The function returns a 'GenomicData' object.
# 	"""
# 	log.info('Entering step "load_dense_matrix"')
# 	source_path = utils.create_chrom_filename(config['process/source_path'], config['process/chrom'])
# 	destination_path = utils.destination_dense_matrices(config['process/destination_path'])
# 	log.debug('- Argument "source_path" filled with "{}"'.format(source_path))
# 	log.debug('- Argument "destination_path" filled with "{}"'.format(destination_path))

# 	try:
# 		vcf = hl.read_matrix_table( sourcePath )
# 		x = [y.get('s') for y in vcf.col.collect()]
# 		log.info('Experiments in loaded VCF: {}'.format(len(x)))
# 		log.info('First and last sample: {} // {}'.format(x[ 0 ], x[ len( x ) - 1 ]))
# 		vcf = vcf.transmute_entries(
# 			sample = hl.struct(
# 				sample = vcf.s,
# 				ad = truncateAt(hl,vcf.AD[ 1 ] / hl.sum( vcf.AD ), "2"),
# 				dp = vcf.DP,
# 				gtInt = vcf.GT,
# 				gt = hl.str(vcf.GT),
# 				gq = vcf.GQ
# 			)
# 		)
# 		vcf = vcf.annotate_rows(
# 			ref = vcf.alleles[ 0 ],
# 			alt = vcf.alleles[ 1 ],
# 			pos = vcf.locus.position,
# 			indel = hl.cond(
# 				(hl.len(vcf.alleles[ 0 ]) != (hl.len(vcf.alleles[ 1 ]))) | (hl.len(vcf.alleles[ 0 ]) != 1) | (hl.len(vcf.alleles[ 0 ]) != 1), True, False 
# 			),
# 			samples_germline = hl.filter(
# 				lambda x: (x.dp > MIN_DP) & (x.gq > MIN_GQ), hl.agg.collect(vcf.sample)
# 			)
# 		)
# 		vcf = vcf.filter_rows(hl.agg.any((vcf.sample.gtInt.is_non_ref()) & (vcf.sample.dp > 10) & (vcf.sample.gq > 20)))
# 		vcf = vcf.key_rows_by(vcf.locus, vcf.alleles).distinct_by_row()
# 		if destination_path is not None:
# 			log.debug('Output VCF file will be saved to "{}". Contents will be overwritten.'.format(destination_path))
# 			vcf = vcf.write(destination_path, overwrite = True)
# 		var = GenomicData()
# 		var.vcf = vcf
# 		var.state = ['dense_matrix']
# 		var.file = [destination_path]
# 		return var
# 	except Exception as ex:
# 		log.debug('Unexpected error during the load of dense matrix "{}"'.format(source_path))
# 		log.error('Unexpected error:\n{}'.format(str(ex)))
# 		log.debug('Stack: {}'.format(str(format_exc())))
# 		sys.exit(2)

#def load_germline(log, hl, nPartitions, source_path, destination_path = None):
def germline(config = None, hl = None, log = None):
	"""Function to load a VCF from POSIX and stores it into HDFS.

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

	Returns
	-------
	The function returns a 'GenomicData' object with the loaded data.
	"""
	self = GenomicData()
	self.state = []
	self.file = []
	if log is None:
		self.log = VoidLog()
	else:
		self.log = log

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering loading step "germline"')

	if config is None:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if hl is None:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(config['process/source_file'], config['process/chrom'])
	source_path = utils.create_chrom_filename(config['process/source_path'], config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_file = utils.create_chrom_filename(config['process/destination_file'], config['process/chrom'])
	destination_path = config['process/destination_path']
	autosave = config['process/autosave']

	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')


	self.data = hl.split_multi_hts(hl.import_vcf(str(source_path), array_elements_required = False, force_bgz = True, min_partitions = config['resources/number_of_partitions']))
	x = [y.get('s') for y in self.data.col.collect()]
	self.log.info('> 	. Experiments in loaded VCF: {}'.format(len(x)))
	self.log.info('> 	. First and last sample: {} // {}'.format(x[ 0 ], x[len(x) - 1]))
	self.data = self.data.transmute_entries(sample = hl.struct(
		sample=self.data.s,
		ad = truncateAt(hl,self.data.AD[ 1 ]/hl.sum(self.data.AD), '2'),
		dp = self.data.DP,
		gtInt = self.data.GT,
		gt = hl.str(self.data.GT),
		gq = self.data.GQ
	)).drop('rsid','qual','filters','info')
	self.data = self.data.annotate_rows(
		ref = self.data.alleles[ 0 ],
		alt = self.data.alleles[ 1 ],
		pos = self.data.locus.position,
		indel = hl.cond((hl.len(self.data.alleles[ 0 ]) != (hl.len(self.data.alleles[ 1 ]))) | (hl.len(self.data.alleles[ 0 ]) != 1) | (hl.len(self.data.alleles[ 0 ]) != 1), True, False),
		samples_germline = hl.filter(lambda x: (x.dp > MIN_DP) & (x.gq > MIN_GQ),hl.agg.collect(self.data.sample))
	) 
	self.data = self.data.annotate_rows(
		freqIntGermline = hl.cond((hl.len(self.data.samples_germline) > 0) | (hl.len(hl.filter(lambda x: x.dp > MIN_DP,self.data.samples_germline)) > 0),
		truncateAt(hl,hl.sum(hl.map(lambda x: x.gtInt.unphased_diploid_gt_index(),self.data.samples_germline))/hl.sum(hl.map(lambda x: 2,hl.filter(lambda x: x.dp > MIN_DP,self.data.samples_germline))),"6"), 0.0)
	).drop('sample')
	self.data = self.data.key_rows_by(self.data.locus, self.data.alleles )

	# if (originPath != ""):
	# 	print ("[INFO]:   . Provided origin path '{}' to be loaded and merged.".format(originPath))
	# 	somatic = hl.read_table(originPath)
	# 	self.data = merge(hl, self.data, somatic)
	# print ("[INFO]:   . Output VCF file will be saved to '{}'".format(destinationPath))
	
	self.state = ['germline'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_germline(destination_path, destination_file)
		self.data.distinct_by_row().write(filename, overwrite = True )
		self.file = [destination_path] + self.file
	return self


# def importSomaticFile(hl, file_path, num_partitions):
#     """ Imports a single somatic vcf file
#         :param HailContext hl: The Hail context
#         :param String file_path: Path from which to import the file
#         :param Int num_partitions: Number of partitions when importing the file
#     """
#     dataset = hl.split_multi_hts(hl.import_vcf(file_path,force_bgz=True,min_partitions=num_partitions)) 
#     return annotateSomatic(hl,dataset)
    
# def importSomatic(hl, originPath, file_paths, destination_path, num_partitions):
#     """ Imports a set of somatic files and merges them into a table. It also merges them with 
#         previously imported germline samples, if any.
#         :param HailContext hl: The hail context
#         :param String originPath: Origin path from which to import previously imported germline samples, if any
#         :param String file_paths: List of file paths from which to import the files
#         :param String destination_path: Path where the loaded variants will be stored
#         :param Int num_partitions: Number of partitions when importing the file
#     """
#     print('[INFO]: Starting process "importSomatic"')
#     nFiles = len(file_paths)
#     print('[INFO]:   . Total number of files to process: {0}'.format(nFiles))
#     print('[INFO]:   . First and last file: {0} / {1}'.format(file_paths[0], file_paths[nFiles - 1]))
#     if(nFiles > 0) :
#         try:
#             # print(file_paths)
#             # print(len(file_paths))
#             tables = [None] * len(file_paths)
#             iteration = 0
#             if (len(tables) == 1):
#                 tables[0] = importSomaticFile(hl,file_paths[0],num_partitions)
#             else:
#                 while (len(tables) > 1):
#                     tmp = []
#                     iteration += 1
#                     print("Iteration ----> " + str(iteration))
#                     for i in range(0, len(tables), 2):
#                         iNext = i+1
#                         if (iteration > 1): 
#                             if (iNext < len(tables)):
#                                 tmp.append(mergeSomatic(hl,tables[i],tables[i+1]))
#                             else:
#                                 tmp.append(tables[i])
#                         else:
#                             table = importSomaticFile(hl,file_paths[i],num_partitions)
#                             if (iNext < len(tables)):
#                                 tableNext = importSomaticFile(hl,file_paths[i+1],num_partitions)
#                                 tmp.append(mergeSomatic(hl,table,tableNext))
#                             else:
#                                 tmp.append(table)
#                     tables = tmp
#             merged = tables[0]
#             if (originPath != ""):
#                 try:
#                     germline = hl.read_table(originPath)
#                     merged = merge(hl,germline,merged)
#                 except Exception:
#                     print('[ERR]: An error was encountered when loading and merging origin content. Was it from germline?')
#                     raise Exception
#             merged.write(destination_path,overwrite=True)
#         except ValueError:
#             print("Error in loading vcf")
#     else:
#         print("Empty file list")
#         if (originPath != ""):
#             germline = hl.read_table(originPath)
#             germline.write(destinationPath,overwrite=True)

# def mergeSomatic(hl, tdataset, tother):
#     """ Merges somatic variants into a single dataset
#         :param HailContext hl:
#         :param Table tdataset: Table with the previosuly merged variants
#         :param Table other: Table to merge with the previously merged tables
#     """
#     joined = tdataset.join(tother,"outer")
#     return joined.transmute(
#         samples_somatic = hl.cond(hl.is_defined(joined.samples_somatic) & hl.is_defined(joined.samples_somatic_1),
#                 joined.samples_somatic.extend(joined.samples_somatic_1),
#                 hl.or_else(joined.samples_somatic,joined.samples_somatic_1)
#                ),
#         was_split = hl.or_else(joined.was_split,joined.was_split_1),
#         a_index = hl.or_else(joined.a_index,joined.a_index_1),
#         ref = hl.or_else(joined.ref,joined.ref_1),
#         alt = hl.or_else(joined.alt,joined.alt_1),
#         pos = hl.or_else(joined.pos,joined.pos_1),
#         indel = hl.or_else(joined.indel,joined.indel_1),
#         freqIntSomatic = hl.or_else(joined.freqIntSomatic,joined.freqIntSomatic_1)
#     )

# def merge(hl, tgermline, tsomatic):
#     joined = tgermline.join(tsomatic,"outer")
#     return joined.transmute(
#         was_split = hl.or_else(joined.was_split,joined.was_split_1),
#         a_index = hl.or_else(joined.a_index,joined.a_index_1),
#         ref = hl.or_else(joined.ref,joined.ref_1),
#         alt = hl.or_else(joined.alt,joined.alt_1),
#         pos = hl.or_else(joined.pos,joined.pos_1),
#         indel = hl.or_else(joined.indel,joined.indel_1)
#     )

# def annotateChrom(hl,chrom):
#     """ Converts input string chromosomes to their integer representation
#         :param Hailcontext hl: The Hail context
#         :param String chrom: String representation of a single chromosome
#     """
#     print("chrom to int function")
#     print(chrom)
#     return (hl.case()
#                     .when(chrom == "MT", "23")
#                     .when(chrom == "X", "24")
#                     .when(chrom == "Y", "25")
#                     .when(chrom == "All", "")
#                     .default(chrom))

# def loadCNV(hl, sourcePath, destinationPath, nPartitions):
#     """ Load CNV data from a tabular formatted input
#         :param HailContext hl: The Hail context
#         :param String sourcePath: The source path from which to import the CNVs
#         :param destinationPath: Path where the loaded variants will be stored
#         :param nPartitions: Number of partitions when importing the file
#     """
#     table = hl.import_table(sourcePath,min_partitions=nPartitions) 
#     table = table.select(
#         table.sample_id,
#         table.start,
#         table.end,
#         table.type,
#         table.cnt,
#         table.chrom,
#         table.bf,
#         table.DGV_goldstd_overlap,
#         table.DGV_goldstd_coordinates,
#         table.genes,
#         table.omim_number,
#         table.omim_phenotype,
#         table.reads_expected,
#         table.reads_observed,
#         table.reads_ratio
#     ) 
#     table.annotate(chrom=annotateChrom(hl,table.chrom),
#                    genes=table.genes.split(",").map(lambda x: hl.struct(gene_name=x)),
#                    intFreqCNV=truncateAt(hl,hl.float(hl.int(table.cnt)/SAMPLES_CNV),"6"),
#                    length=hl.abs(hl.int(table.end)-hl.int(table.start))) \
#          .write(destinationPath,overwrite=True) 