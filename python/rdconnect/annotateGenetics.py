
import os
import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData
from rdconnect.classException import *
from rdconnect.utils import _check_class_and_config

def _transcript_annotations(hl, annotations):
	""" Transcript level annotations for VEP 
		:param Hailcontext hl: The Hail context
		:param HailTable: Previously annotated data (variant level)
	"""
	return hl.map(lambda x: 
		hl.struct(
			gene_name = x.gene_symbol,
			effect_impact = x.impact,
			transcript_id = x.transcript_id,
			effect = hl.delimit(x.consequence_terms, ','),
			gene_id = x.gene_id,
			functional_class = 'transcript',
			amino_acid_length = '',
			codon_change = x.hgvsc.replace('.*:', ''),
			amino_acid_change = x.hgvsp.replace('.*:', ''),
			exon_rank = x.exon,
			transcript_biotype = x.biotype,
			gene_coding = hl.str(x.cds_start)
		),
	annotations)


def _intergenic_annotations(hl, annotations):
	""" Transcript level annotations for VEP 
		:param Hailcontext hl: The Hail context
		:param HailTable: Previously annotated data (variant level)
	"""
	return hl.map(lambda x: 
		hl.struct(
			gene_name = '',
			effect_impact = x.impact,
			transcript_id = '',
			effect = hl.delimit(x.consequence_terms, ','),
			gene_id = '',
			functional_class = 'intergenic_region',
			amino_acid_length = '0',
			codon_change = '',
			amino_acid_change = '',
			exon_rank = '',
			transcript_biotype = '',
			gene_coding = ''
		),
	annotations)

def vep(self = None, config = None, hl = None, log = None):
	"""Annotates given genetic dataset with VEP annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
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
	The function returns a 'GenomicData' object annotated with VEP.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "VEP"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_path = self.config['process/destination_path']
	vep_config = self.config['annotation/clean/vep']
	autosave = self.config['process/autosave']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "vep_config" filled with "{}"'.format(vep_config))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')

	if 'data' not in vars(self):
		self.log.info('Loading genomic data from "source_path"')
		self.data = self.hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = self.hl.vep(self.data, vep_config)
	self.data = self.data.annotate_rows(
		effs = self.hl.cond(
			self.hl.is_defined(self.data.vep.transcript_consequences),
			_transcript_annotations(self.hl,self.data.vep.transcript_consequences),
			_intergenic_annotations(self.hl,self.data.vep.intergenic_consequences)
		),
		rs = self.data.vep.colocated_variants[0].id
	)
	self.data = self.data.drop('vep')

	self.state = ['VEP'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_vep(destination_path, source_file)
		self.data.write(filename, overwrite = True)
		self.file = [destination_path] + self.file
	return self


def _mt_pred_annotations(hl, annotations):
	""" Annotations for dbNSFP
		:param Hailcontext hl: The Hail context
		:param HailTable: Previously annotated data (variant level)
	"""
	arr = annotations.MutationTaster_pred.split(";")
	return (hl.case()
		.when(arr.contains("A"),"A")
		.when(arr.contains("D"),"D")
		.when(arr.contains("N"),"N")
		.default(""))


def _polyphen_pred_annotations(hl, annotations):
	""" Annotations for dbNSFP
		:param Hailcontext hl: The Hail context
		:param HailTable: Previously annotated data (variant level)
	"""
	arr = annotations.Polyphen2_HDIV_pred.split(";")
	return (hl.case()
		.when(arr.contains("D"),"D")
		.when(arr.contains("P"),"P")
		.when(arr.contains("B"),"B")
		.default("")
	)


def _sift_pred_annotations(hl, annotations):
	""" Annotations for dbNSFP
		:param Hailcontext hl: The Hail context
		:param HailTable: Previously annotated data (variant level)
	"""
	arr = annotations.SIFT_pred
	return (hl.case()
		.when(arr.contains("D"),"D")
		.when(arr.contains("T"),"T")
		.default("")
	)


def _truncateAt(hl, n, p):
	""" Formats a input number to 'p' decimals
		:param Hailcontext hl: The Hail context
		:param String n: Number to format
		:param String p: Decimal precision
	"""
	return hl.float(hl.int((10 ** hl.int(p) * n))) / (10 ** hl.int(p))


def _removeDot(hl, n, precision):
	""" Formats an input number to 'p' decimals, or sets it to 0 if it's a dot annotation
		:param HailContext h: The Hail context
		:param String n: Number to format
		:param String p: Decimal precision
	"""
	return hl.cond(n.startswith('.'),0.0,_truncateAt(hl,hl.float(n),precision))


def dbnsfp(self = None, config, hl = None, log = None):
	"""Annotates given genetic dataset with dbSNFP annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
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
	The function returns a 'GenomicData' object annotated with dbSNFP.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "dbSNFP"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_path = self.config['process/destination_path']
	dbnsfp_path = utils.create_chrom_filename(self.config['annotation/clean/dbNSFP'], self.config['process/chrom'])
	autosave = self.config['process/autosave']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "dbnsfp_path" filled with "{}"'.format(dbnsfp_path))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')

	dbnsfp = self.hl.read_table(dbnsfp_path)
	if 'data' not in vars(self):
		self.log.info('Loading genomic data from "source_path"')
		self.data = self.hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = self.data.annotate_rows(
		gp1_asn_af = self.hl.or_else(_removeDot(self.hl, dbnsfp[self.data.locus, self.data.alleles].Gp1_ASN_AF1000, '6'), 0.0),
		gp1_eur_af = self.hl.or_else(_removeDot(self.hl, dbnsfp[self.data.locus, self.data.alleles].Gp1_EUR_AF1000, '6'), 0.0),
		gp1_afr_af = self.hl.or_else(_removeDot(self.hl, dbnsfp[self.data.locus, self.data.alleles].Gp1_AFR_AF1000, '6'), 0.0),
		gp1_af = self.hl.or_else(_removeDot(self.hl, dbnsfp[self.data.locus, self.data.alleles].Gp1_AF1000, '6'), 0.0),
		gerp_rs = dbnsfp[self.data.locus, self.data.alleles].GERP_RS,
		mt = self.hl.or_else(self.hl.max(dbnsfp[self.data.locus, self.data.alleles].MutationTaster_score.split(';').map(lambda x:_removeDot(self.hl, x, '4'))), 0.0),
		mutationtaster_pred = _mt_pred_annotations(self.hl, dbnsfp[self.data.locus, self.data.alleles]),
		phyloP46way_placental = _removeDot(self.hl, dbnsfp[self.data.locus, self.data.alleles].phyloP46way_placental, '4'),
		polyphen2_hvar_pred = _polyphen_pred_annotations(self.hl, dbnsfp[self.data.locus, self.data.alleles]),
		polyphen2_hvar_score = self.hl.or_else(self.hl.max(dbnsfp[self.data.locus, self.data.alleles].Polyphen2_HVAR_score.split(';').map(lambda x: _removeDot(self.hl, x, '4'))), 0.0),
		sift_pred = _sift_pred_annotations(self.hl, dbnsfp[self.data.locus, self.data.alleles]),
		sift_score = self.hl.or_else(self.hl.max(dbnsfp[self.data.locus, self.data.alleles].SIFT_score.split(';').map(lambda x: _removeDot(self.hl, x, '4'))), 0.0),
		cosmic_id = dbnsfp[self.data.locus, self.data.alleles].COSMIC_ID)

	self.state = ['dbNSFP'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_dbnsfp(destination_path, source_file)
		self.data.write(filename, overwrite = True)
		self.file = [destination_path] + self.file
	return self


def cadd(self = None, config = None, hl = None, log = None):
	"""Annotates given genetic dataset with CADD annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
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
	The function returns a 'GenomicData' object annotated with CADD.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "CADD"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_path = self.config['process/destination_path']
	cad_path = utils.create_chrom_filename(self.config['annotation/clean/cadd'], self.config['process/chrom'])
	autosave = self.config['process/autosave']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "cad_path" filled with "{}"'.format(cad_path))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')

	cadd = self.hl.split_multi_hts(self.hl.read_matrix_table(cad_path)) \
		.rows() \
		.key_by('locus', 'alleles')

	if 'data' not in vars(self):
		self.log.info('Loading genomic data from "source_path"')
		self.data = self.hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = self.data.annotate_rows(cadd_phred = cadd[self.data.locus, self.data.alleles].info.CADD13_PHRED[cadd[self.data.locus, self.data.alleles].a_index-1])

	self.state = ['CADD'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_cadd(destination_path, source_file)
		self.data.write(filename, overwrite = True)
		self.file = [destination_path] + self.file
	return self


def _clinvar_filtering(hl, annotation, is_filter_field):
	""" Returns the ClinVar annotations to apply
		:param HailContext hl: The Hail context
		:param Annotation annotation: Annotations to apply
		:param Boolean is_filter_field: Whether the annotations is for ClinVar filtering or informative
	"""
	clin_sigs = hl.dict([
		('Uncertain_significance', 'VUS'),
		('not_provided', 'NA'),
		('Benign', 'B'),
		('Likely_benign', 'LB'),
		('Likely_pathogenic', 'LP'),
		('Pathogenic', 'P'),
		('drug_response', 'Drug'),
		('histocompatibility', 'Histo'),
		('Conflicting_interpretations_of_pathogenicity', 'C'),
		('Affects', 'Other'),
		('risk_factor', 'Other'),
		('association', 'Other'),
		('protective', 'Other'),
		('other', 'Other')
	])
	filtered = None
	if is_filter_field:
		filtered = hl.map(lambda z: hl.cond(clin_sigs.contains(z), hl.struct(clnsig=clin_sigs[z]), hl.struct(clnsig="-1")), annotation)
		filtered = hl.filter(lambda e: e['clnsig'] != '-1', filtered)    
	else: 
		filtered = hl.map(lambda z: hl.cond(clin_sigs.contains(z), clin_sigs[z], '-1'), annotation)
		filtered = hl.filter(lambda e: e != '-1', filtered)  
	return filtered


def _clinvar_preprocess(hl, annotation, is_filter_field):
	""" Preprocesses a Clinvar annotation expression
		:param Hailcontext hl: The Hail context
		:param Annotation annotation: Annotations to apply
		:param Boolean is_filter_field: Whether the annotations is for Clinvar filtering or informative
	"""
	preprocessed = hl.flatmap(lambda x: x.replace('\\/',',')
		.replace('\\:',',') \
		.replace('\\|',',') \
		.split(','), annotation)
	preprocessed = hl.map(lambda y: hl.cond(y[0] == '_', y[1:], y), preprocessed)
	return _clinvar_filtering(hl,preprocessed,is_filter_field)


def clinvar(self = None, config = None, hl = None, log = None):
	"""Annotates given genetic dataset with ClinVar annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
	config: ConfigFile, optional
		Configuration for this step of the pipeline.
	hl: context, optional
		HAIL context.
	log: logger, optional
		A logger to have track of the steps used in the loading process.

	Returns
	-------
	The function returns a 'GenomicData' object annotated with ClinVar.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "ClinVar"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_path = self.config['process/destination_path']
	clinvar_path = utils.create_chrom_filename(self.config['annotation/clean/clinvar'], self.config['process/chrom'])
	autosave = self.config['process/autosave']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "clinvar_path" filled with "{}"'.format(clinvar_path))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')

	clinvar = self.hl.split_multi_hts(self.hl.read_matrix_table(clinvar_path)) \
		.rows() \
		.key_by('locus', 'alleles')

	if 'data' not in vars(self):
		self.log.info('Loading genomic data from "source_path"')
		self.data = self.hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = self.data.annotate_rows(
		clinvar_id = self.hl.cond(self.hl.is_defined(clinvar[self.data.locus, self.data.alleles].info.CLNSIG[clinvar[self.data.locus, self.data.alleles].a_index-1]), clinvar[self.data.locus, self.data.alleles].rsid, clinvar[self.data.locus, self.data.alleles].info.CLNSIGINCL[0].split(':')[0]),
		clinvar_clnsigconf = self.hl.delimit(clinvar[self.data.locus, self.data.alleles].info.CLNSIGCONF),
		clinvar_clnsig = self.hl.cond(self.hl.is_defined(clinvar[self.data.locus, self.data.alleles].info.CLNSIG[clinvar[self.data.locus, self.data.alleles].a_index-1]), self.hl.delimit(_clinvar_preprocess(self.hl,clinvar[self.data.locus, self.data.alleles].info.CLNSIG,False), "|"), self.hl.delimit(_clinvar_preprocess(self.hl,clinvar[self.data.locus, self.data.alleles].info.CLNSIGINCL, False), "|")),
		clinvar_filter = self.hl.cond(self.hl.is_defined(clinvar[self.data.locus, self.data.alleles].info.CLNSIG[clinvar[self.data.locus, self.data.alleles].a_index-1]), _clinvar_preprocess(self.hl,clinvar[self.data.locus, self.data.alleles].info.CLNSIG,True), _clinvar_preprocess(self.hl, clinvar[self.data.locus, self.data.alleles].info.CLNSIGINCL, True))
	)

	self.state = ['ClinVar'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_clinvar(destination_path, source_file)
		self.data.write(filename, overwrite = True)
		self.file = [destination_path] + self.file
	return self



def gnomADEx(self = None, config = None, hl = None, log = None):
	"""Annotates given genetic dataset with gnomeAD exome annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
	config: ConfigFile, optional
		Configuration for this step of the pipeline.
	hl: context, optional
		HAIL context.
	log: logger, optional
		A logger to have track of the steps used in the loading process.

	Returns
	-------
	The function returns a 'GenomicData' object annotated with ClinVar.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = _check_class_and_config(self, config, hl, log)
	self.log.info('Entering annotation step "gnomeAD Exome"')

	if not isConfig:
		self.log.error('No configuration was provided')
		raise NoConfigurationException('No configuration was provided')

	if not isHl:
		self.log.error('No pointer to HAIL module was provided')
		raise NoHailContextException('No pointer to HAIL module was provided')

	source_file = utils.create_chrom_filename(self.config['process/source_file'], self.config['process/chrom'])
	source_path = utils.create_chrom_filename(self.config['process/source_path'], self.config['process/chrom'])
	source_path = os.path.join(source_path, source_file)
	destination_path = self.config['process/destination_path']
	gnomeAdEx_path = utils.create_chrom_filename(self.config['annotation/clean/exomesGnomad'], self.config['process/chrom'])
	autosave = self.config['process/autosave']

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
	self.log.debug('> Argument "gnomeAdEx_path" filled with "{}"'.format(gnomeAdEx_path))
	self.log.debug('> Argument "autosave" was set' if autosave else '> Argument "autosave" was not set')

	gnomad = self.hl.split_multi_hts(self.hl.read_matrix_table(gnomeAdEx_path)) \
		.rows() \
		.key_by('locus','alleles')

	self.data = self.data.annotate_rows(
		gnomad_af = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AF[gnomad[self.data.locus, self.data.alleles].a_index-1]), gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AF[gnomad[self.data.locus, self.data.alleles].a_index-1], 0.0),
		gnomad_ac = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AC[gnomad[self.data.locus, self.data.alleles].a_index-1]), gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AC[gnomad[self.data.locus, self.data.alleles].a_index-1], 0.0),
		gnomad_an = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AN),gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AN, 0.0),
		gnomad_af_popmax = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AF_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1]), gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AF_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1], 0.0),
		gnomad_ac_popmax = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AC_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1]), gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AC_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1], 0.0),
		gnomad_an_popmax = self.hl.cond(self.hl.is_defined(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AN_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1]), gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_AN_POPMAX[gnomad[self.data.locus, self.data.alleles].a_index-1], 0.0),
		gnomad_filter = self.hl.cond(gnomad[self.data.locus, self.data.alleles].info.gnomAD_Ex_filterStats == 'Pass','PASS','non-PASS')
	)

	self.state = ['gnomeADEx'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_clinvar(destination_path, source_file)
		self.data.write(filename, overwrite = True)
		self.file = [destination_path] + self.file
	return self
