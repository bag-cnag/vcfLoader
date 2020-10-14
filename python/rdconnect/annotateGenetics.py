
import utils
from classGenome import GenomicData


def VEP(log, hl, var, source_path, destination_path, vep_config, autosave):
	"""Annotates given genetic dataset with VEP annotations.

	Parameters
	----------
	log: logger, mandatory
		A logger to have track of the steps used in the loading process.
	hl: context, mandatory
		HAIL context.
	var: GenomicData, mandatory
		Set it to None to load the dataset from 'source_path'. If a GenomicData
		is assigned to this argument, no set is loaded from 'source_path' and
		the argument is ignored.
	source_path: str, mandatory
		String used to load the 'MatrixTable' from HAIL set to be annotated.
	destination_path: str, mandatory
		String used to save the annotated 'MatrixTable'. Only used of 
		'autosave' is set to True, ignored otherwise.
	vep_config: str, mandatory
		VEP configuration path
	autosave: bool, mandatory
		If set to True, the annotated 'MatrixTable' is saved as file using the
		template in 'destination_path'.
	Returns
	-------
	The function returns a 'GenomicData' object.
	"""
	log.info('Entering annotation step "VEP"')
	if var is None:
		log.debug('- Argument "var" was not set')
	else:
		log.debug('- Argument "var" was set')
	log.debug('- Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('- Argument "destination_path" filled with "{}"'.format(destination_path))
	log.debug('- Argument "vep_config" filled with "{}"'.format(vep_config))
	if autosave:
		log.debug('- Argument "autosave" was set')
	else:
		log.debug('- Argument "autosave" was not set')

	if var is None:
		var = GenomicData()
		var.data = hl.methods.read_matrix_table(source_path)
		var.state = []
		var.file = []

	var.data = hl.vep(var.data, vep_config)
	var.data = var.data.annotate_rows(
		effs = hl.cond(
			hl.is_defined(var.data.vep.transcript_consequences),
			transcript_annotations(hl,var.data.vep.transcript_consequences),
			intergenic_annotations(hl,var.data.vep.intergenic_consequences)
		),
		rs = var.data.vep.colocated_variants[0].id
	)
	var.data = var.data.drop("vep")

	var.state = ['VEP'] + var.state
	if autosave:
		filename = utils.destination_vep(destination_path)
		var.data.write(destination_path, overwrite = True)
		var.file = [destination_path] + var.file
	return var

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


def dbNSFP(log, hl, var, source_path, destination_path, dbnsfp_path, autosave):
	""" Adds dbNSFP annotations to variants.
		:param HailContext hl: The Hail context
		:param VariantDataset variants: The variants to annotate
		:param string dbnsfpPath: Path were the dbNSFP table can be found
		:param string destinationPath: Path were the new annotated dataset can be found
	"""
	log.info('Entering annotation step "VEP"')
	if var is None:
		log.debug('- Argument "var" was not set')
	else:
		log.debug('- Argument "var" was set')
	log.debug('- Argument "source_path" filled with "{}"'.format(source_path))
	log.debug('- Argument "destination_path" filled with "{}"'.format(destination_path))
	log.debug('- Argument "dbnsfp_path" filled with "{}"'.format(dbnsfp_path))
	if autosave:
		log.debug('- Argument "autosave" was set')
	else:
		log.debug('- Argument "autosave" was not set')

	dbnsfp = hl.read_table(dbnsfp_path)
	if var is None:
		var = GenomicData()
		var.data = hl.methods.read_matrix_table(source_path)
		var.state = []
		var.file = []

	var.data = var.data.annotate_rows(
		gp1_asn_af=hl.or_else(_removeDot(hl,dbnsfp[var.data.locus, var.data.alleles].Gp1_ASN_AF1000,"6"), 0.0),
		gp1_eur_af=hl.or_else(_removeDot(hl,dbnsfp[var.data.locus, var.data.alleles].Gp1_EUR_AF1000,"6"), 0.0),
		gp1_afr_af=hl.or_else(_removeDot(hl,dbnsfp[var.data.locus, var.data.alleles].Gp1_AFR_AF1000,"6"), 0.0),
		gp1_af=hl.or_else(_removeDot(hl,dbnsfp[var.data.locus, var.data.alleles].Gp1_AF1000,"6"), 0.0),
		gerp_rs=dbnsfp[var.data.locus, var.data.alleles].GERP_RS,
		mt=hl.or_else(hl.max(dbnsfp[var.data.locus, var.data.alleles].MutationTaster_score.split(";").map(lambda x:_removeDot(hl,x,"4"))),0.0),
		mutationtaster_pred=_mt_pred_annotations(hl,dbnsfp[var.data.locus, var.data.alleles]),
		phyloP46way_placental=_removeDot(hl,dbnsfp[var.data.locus, var.data.alleles].phyloP46way_placental,"4"),
		polyphen2_hvar_pred=_polyphen_pred_annotations(hl,dbnsfp[var.data.locus, var.data.alleles]),
		polyphen2_hvar_score=hl.or_else(hl.max(dbnsfp[var.data.locus, var.data.alleles].Polyphen2_HVAR_score.split(";").map(lambda x: _removeDot(hl,x,"4"))),0.0),
		sift_pred=_sift_pred_annotations(hl,dbnsfp[var.data.locus, var.data.alleles]),
		sift_score=hl.or_else(hl.max(dbnsfp[var.data.locus, var.data.alleles].SIFT_score.split(";").map(lambda x: _removeDot(hl,x,"4"))),0.0),
		cosmic_id=dbnsfp[var.data.locus, var.data.alleles].COSMIC_ID)

	var.state = ['dbNSFP'] + var.state
	if autosave:
		filename = utils.destination_dbnsfp(destination_path)
		var.data.write(destination_path, overwrite = True)
        var.file = [destination_path] + var.file
	return var