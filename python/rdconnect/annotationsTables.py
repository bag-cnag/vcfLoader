
import os
import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData
from rdconnect.classException import *
from rdconnect.utils import check_class_and_config


def dbsnfp(config, hl, log):
	"""Loads into the HDFS in HAIL compatible format dbSNFP annotation.

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
	self, isConfig, isHl = check_class_and_config(None, config, hl, log)
	self.log.info('Entering importing step "dbSNFP"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['annotation/raw/dbNSFP'], self.config['process/chrom'])
	destination_file = utils.create_chrom_filename(self.config['annotation/clean/dbNSFP'], self.config['process/chrom'])

	self.log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))

	table = hl.import_table(sourcePath, min_partitions = self.config['resources/number_of_partitions']) \
		.rename({
			'#chr': 'chr',
			'pos(1-coor)': 'pos',
			'1000Gp1_AF':'Gp1_AF1000',
			'1000Gp1_AC':'Gp1_AC1000',
			'1000Gp1_EUR_AF':'Gp1_EUR_AF1000',
			'1000Gp1_ASN_AF':'Gp1_ASN_AF1000',
			'1000Gp1_AFR_AF':'Gp1_AFR_AF1000',
			'ESP6500_EA_AF ':'ESP6500_EA_AF',
			'GERP++_RS':'GERP_RS'
		}
	)
	table = table.annotate(locus = hl.locus(table.chr,hl.int(table.pos)), alleles = [table.ref, table.alt]) 
	table = table.select(
		table.locus,
		table.alleles,
		table.Gp1_AF1000,
		table.Gp1_EUR_AF1000,
		table.Gp1_ASN_AF1000,
		table.Gp1_AFR_AF1000,
		table.GERP_RS,
		table.MutationTaster_score,
		table.MutationTaster_pred,
		table.phyloP46way_placental,
		table.Polyphen2_HDIV_pred,
		table.Polyphen2_HVAR_score,
		table.SIFT_pred,
		table.SIFT_score,
		table.COSMIC_ID
	)
	table.key_by('locus','alleles') \
		.write(destination_file, overwrite = True) 

def cgi(config, hl, log):
	"""Loads into the HDFS in HAIL compatible format CGI annotation.

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
	self, isConfig, isHl = check_class_and_config(None, config, hl, log)
	self.log.info('Entering importing step "CGI"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['annotation/raw/CGI'], self.config['process/chrom'])
	destination_file = utils.create_chrom_filename(self.config['annotation/clean/CGI'], self.config['process/chrom'])

	self.log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))

	table = hl.import_table(source_file, min_partitions = self.config['resources/number_of_partitions']) \
		.rename({'#CHRM': 'chr'})
	table.annotate(locus = hl.locus(table.chr, hl.int(table.POS)), alleles = [table.REF,table.ALT]) \
		.key_by('locus', 'alleles') \
		.write(destination_file, overwrite = True) 

def cadd(config, hl, log):
	"""Loads into the HDFS in HAIL compatible format CADD annotation.

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
	self, isConfig, isHl = check_class_and_config(None, config, hl, log)
	self.log.info('Entering importing step "CADD"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['annotation/raw/cadd'], self.config['process/chrom'])
	destination_file = utils.create_chrom_filename(self.config['annotation/clean/cadd'], self.config['process/chrom'])

	self.log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))

	vcf(self.hl, source_file, destination_file, self.config['resources/number_of_partitions'])


def clinvar(config, hl, log):
	"""Loads into the HDFS in HAIL compatible format ClinVar annotation.

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
	self, isConfig, isHl = check_class_and_config(None, config, hl, log)
	self.log.info('Entering importing step "ClinVar"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['annotation/raw/clinvar'], self.config['process/chrom'])
	destination_file = utils.create_chrom_filename(self.config['annotation/clean/clinvar'], self.config['process/chrom'])

	self.log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))

	vcf(self.hl, source_file, destination_file, self.config['resources/number_of_partitions'])


def gnomADEx(config, hl, log):
	"""Loads into the HDFS in HAIL compatible format gnomeAD Exon annotation.

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
	self, isConfig, isHl = check_class_and_config(None, config, hl, log)
	self.log.info('Entering importing step "ClinVar"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['annotation/raw/exomesGnomad'], self.config['process/chrom'])
	destination_file = utils.create_chrom_filename(self.config['annotation/clean/exomesGnomad'], self.config['process/chrom'])

	self.log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
	self.log.debug('> Argument "destination_file" filled with "{}"'.format(destination_file))

	vcf(self.hl, source_file, destination_file, self.config['resources/number_of_partitions'])


def vcf(hl, sourcePath, destinationPath, nPartitions):
	""" Imports annotations vcfs
		:param HailContext hl: The Hail context
		:param String sourcePath: Annotation vcf path
		:param String destinationPath: Path where the loaded annotation file will be put
		:param String nPartitions: Number of partitions
	"""
	hl.import_vcf(sourcePath,min_partitions=nPartitions,skip_invalid_loci=True) \
		.write(destinationPath,overwrite=True)
