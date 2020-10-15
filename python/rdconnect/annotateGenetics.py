
import rdconnect.utils as utils
from rdconnect.classGenome import GenomicData


def vep(self, config, hl, log = None):
	"""Annotates given genetic dataset with VEP annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
	config: ConfigFile, mandatory
		Configuration for this step of the pipeline.
	hl: context, mandatory
		HAIL context.
	log: logger, optional
		A logger to have track of the steps used in the loading process.

	Returns
	-------
	The function returns a 'GenomicData' object.
	"""
	if log is not None: 
		log.info('Entering annotation step "VEP"')

	source_file = utils.create_chrom_filename(config['process/source_file'], config['process/chrom'])
	source_path = utils.create_chrom_filename(config['process/source_path'], config['process/chrom'])
	destination_path = config['process/destination_path']
	vep_config = config['annotation/clean/vep_config']
	autosave = config['process/autosave']
	print('autosave', autosave)

	if self is None and log is not None:
		log.debug('> Argument "self" was not set')
	if self is not None and log is not None:
		log.debug('> Argument "self" was set')
	if log is not None: 
		log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
		log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
		log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
		log.debug('> Argument "vep_config" filled with "{}"'.format(vep_config))
	if autosave and log is not None:
		log.debug('> Argument "autosave" was set')
	if not autosave and log is not None:
		log.debug('> Argument "autosave" was not set')

	if self is None:
		self = GenomicData()
		self.data = hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = hl.vep(self.data, vep_config)
	self.data = self.data.annotate_rows(
		effs = hl.cond(
			hl.is_defined(self.data.vep.transcript_consequences),
			transcript_annotations(hl,self.data.vep.transcript_consequences),
			intergenic_annotations(hl,self.data.vep.intergenic_consequences)
		),
		rs = self.data.vep.colocated_selfiants[0].id
	)
	self.data = self.data.drop('vep')

	self.state = ['VEP'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_vep(destination_path, source_file)
		self.data.write(destination_path, overwrite = True)
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


def dbnsfp(self, config, hl, log = None):
	"""Annotates given genetic dataset with VEP annotations.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
	config: ConfigFile, mandatory
		Configuration for this step of the pipeline.
	hl: context, mandatory
		HAIL context.
	log: logger, optional
		A logger to have track of the steps used in the loading process.

	Returns
	-------
	The function returns a 'GenomicData' object.
	"""
	if log is not None:
		log.info('Entering annotation step "VEP"')

	source_file = utils.create_chrom_filename(config['process/source_file'], config['process/chrom'])
	source_path = utils.create_chrom_filename(config['process/source_path'], config['process/chrom'])
	destination_path = config['process/destination_path']
	dbnsfp_path = config['annotation/clean/dbNSFP']
	autosave = config['process/autosave']

	if self is None and log is not None:
		log.debug('> Argument "self" was not set')
	if self is not None and log is not None:
		log.debug('> Argument "self" was set')
	if log is not None:
		log.debug('> Argument "source_file" filled with "{}"'.format(source_file))
		log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
		log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))
		log.debug('> Argument "dbnsfp_path" filled with "{}"'.format(dbnsfp_path))
		
	if autosave and log is not None:
		log.debug('> Argument "autosave" was set')
	if not autosave and log is not None:
		log.debug('> Argument "autosave" was not set')

	dbnsfp = hl.read_table(dbnsfp_path)
	if self is None:
		self = GenomicData()
		self.data = hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	self.data = self.data.annotate_rows(
		gp1_asn_af=hl.or_else(_removeDot(hl,dbnsfp[self.data.locus, self.data.alleles].Gp1_ASN_AF1000,"6"), 0.0),
		gp1_eur_af=hl.or_else(_removeDot(hl,dbnsfp[self.data.locus, self.data.alleles].Gp1_EUR_AF1000,"6"), 0.0),
		gp1_afr_af=hl.or_else(_removeDot(hl,dbnsfp[self.data.locus, self.data.alleles].Gp1_AFR_AF1000,"6"), 0.0),
		gp1_af=hl.or_else(_removeDot(hl,dbnsfp[self.data.locus, self.data.alleles].Gp1_AF1000,"6"), 0.0),
		gerp_rs=dbnsfp[self.data.locus, self.data.alleles].GERP_RS,
		mt=hl.or_else(hl.max(dbnsfp[self.data.locus, self.data.alleles].MutationTaster_score.split(";").map(lambda x:_removeDot(hl,x,"4"))),0.0),
		mutationtaster_pred=_mt_pred_annotations(hl,dbnsfp[self.data.locus, self.data.alleles]),
		phyloP46way_placental=_removeDot(hl,dbnsfp[self.data.locus, self.data.alleles].phyloP46way_placental,"4"),
		polyphen2_hvar_pred=_polyphen_pred_annotations(hl,dbnsfp[self.data.locus, self.data.alleles]),
		polyphen2_hvar_score=hl.or_else(hl.max(dbnsfp[self.data.locus, self.data.alleles].Polyphen2_HVAR_score.split(";").map(lambda x: _removeDot(hl,x,"4"))),0.0),
		sift_pred=_sift_pred_annotations(hl,dbnsfp[self.data.locus, self.data.alleles]),
		sift_score=hl.or_else(hl.max(dbnsfp[self.data.locus, self.data.alleles].SIFT_score.split(";").map(lambda x: _removeDot(hl,x,"4"))),0.0),
		cosmic_id=dbnsfp[self.data.locus, self.data.alleles].COSMIC_ID)

	self.state = ['dbNSFP'] + self.state
	if autosave and destination_path != '':
		filename = utils.destination_dbnsfp(destination_path, source_file)
		self.data.write(destination_path, overwrite = True)
		self.file = [destination_path] + var.file
	return self