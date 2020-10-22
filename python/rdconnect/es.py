
import os
from rdconnect.classException import *
import rdconnect.utils as utils

def transform(self = None, config = None, hl = None, log = None):
	"""Transforms a given dataset into the dataframe format for ElasticSearch.

	If no 'GenomicData' is given, it is loaded from 'config' using 
	'process/source_path' and 'process/source_file'. The data is written into
	'process/destination_path' and using as a filename's:

		chrom=[chromosome]

	The function returns the same object given in 'self'. If not 'self' is 
	provided and data is loaded from file, the function returns that data
	as 'GenomicData' object before its transformation to ElasticSearch 
	dataframe.


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
	The function returns a 'GenomicData' object the data before it was
	transformed for ElasticSearch.
	"""
	isSelf = True
	if self is None:
		isSelf = False

	self, isConfig, isHl = utils._check_class_and_config(self, config, hl, log)
	self.log.info('Entering transform step')

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

	self.log.debug('> Argument "self" was set' if isSelf else '> Argument "self" was not set')
	self.log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
	self.log.debug('> Argument "destination_path" filled with "{}"'.format(destination_path))

	if 'data' not in vars(self):
		self.log.info('Loading genomic data from "source_path"')
		self.data = hl.methods.read_matrix_table(source_path)
		self.state = []
		self.file = []

	out_file = os.path.join(destination_path + 'chrom=' + str(self.config['process/chrom']))
	vcf = self.data.rows()
	vcf.key_by(vcf.locus, vcf.alleles).distinct()
	vcf.to_spark() \
		.drop('locus.contig', 'locus.position', 'alleles') \
		.write.mode('overwrite').save(out_file)

	if 'flags' not in vars(self):
		self.flags = {'transform': (True, out_file)}

	return self




def _index_exists(host, port, index_name, user, pwd):
    sts = requests.head('http://{}:{}/{}'.format(host, port, index_name), auth=(user, pwd))
    print('http://{}:{}/{}'.format(host, port, index_name))
    print(sts.status_code)
    print(sts.text)
    return sts.status_code == 200

def _create_index(host,port,index_name,data,user,pwd):
	url = "http://" + host + ":" + port + "/" + index_name
	headers = {'Content-Type': 'application/json'}
	response = requests.put(url,data=data,headers=headers,auth=(user,pwd))
	return response.status_code


def create_index_snv(self = None, config = None, log = None):
	"""Creates an index for SNV data into ElasticSearch.

	Location, data access credentials and index name is looked first into the
	'GenomicData' object, If it is not given, it is looked into the 
	'ConfigFile'.

	Parameters
	----------
	self: GenomicData, mandatory
		Set it to None forces the function to used the provided configuration. 
		If set, it is used to get the access to the ElasicSearch.
	config: ConfigFile, optional
		Configuration for this step of the pipeline. If not provided or set to
		None the configuration is looked into the GenomicData in self.
	log: logger, optional
		A logger to have track of the steps used in the loading process. If not
		provided or set to None the logger is looked into the GenomicData in 
		self. If no logger is in the provided nor in the GenomicData, then no
		log is performed.

	Returns
	-------
	If a 'GenomicData' is provided, it returns the genomic data with a set flag
	indicating the results of the index creation:

		(status code is 200?, index' name, status code obtained)

	If no 'GenomicData' is given to the function, it returns the status code
	obtained from the query.
	"""
	if self is None and config is None:
		raise NoConfigurationException('No configuration was provided in form of GenomicData nor ConfigFile (both were None)') 

	if self is not None and 'config' in vars(self) and config is None:
		config = self.config

	if self is not None and not 'config' in vars(self) and config is None:
		raise NoConfigurationException('No configuration was provided in form of ConfigFile (was None) and provided GenomicData does not contain configuration') 

	host = self.config['resources/elasticsearch/host']
	port = self.config['resources/elasticsearch/port']
	index_name = self.config['resources/elasticsearch/index_name']
	user = self.config['resources/elasticsearch/user']
	pwd = self.config['resources/elasticsearch/pwd']
	num_shards = self.config['resources/elasticsearch/num_shards']
	num_replicas = self.config['resources/elasticsearch/num_replicas']
	version = self.config['resources/elasticsearch/type']

	data = """
		{"settings":{"index":{"number_of_shards":""" + num_shards + ""","number_of_replicas":""" + num_replicas + """, "refresh_interval":"-1"}} ,"mappings":{"""+"\"" + version + "\""+""" :{
			"properties":{
				"chrom":{"type":"integer","index":"true"}
				,"pos":{"type":"integer","index":"true"}
				,"ref":{"type":"keyword","index":"false"}
				,"alt":{"type":"keyword","index":"false"}
				,"indel":{"type":"keyword","index":"true"}
				,"freqInt":{"type":"float"}
				,"rs":{"type":"keyword", "index":"true"}
				,"cadd_phred":{"type":"float","index":"true"}
				,"gerp_rs":{"type":"keyword","index":"false"}
				,"mt":{"type":"float","index":"false"}
				,"mutationtaster_pred":{"type":"keyword"}
				,"phylop46way_placental":{"type":"keyword","index":"false"}
				,"polyphen2_hvar_pred":{"type":"keyword"}
				,"polyphen2_hvar_score":{"type":"float","index":"false"}
				,"sift_pred":{"type":"keyword"}
				,"sift_score":{"type":"float","index":"false"}
				,"siphy_29way_pi":{"type":"keyword","index":"false"}
				,"UMD":{"type":"keyword"}
				,"clinvar_clnsig":{"type":"keyword","index":"false"}
				,"clinvar_clnsigconf":{"type":"keyword","index":"false"}
				,"clinvar_id":{"type":"integer","index":"false"}
				,"gp1_afr_af":{"type":"float","index":"false"}
				,"gp1_asn_af":{"type":"float","index":"false"}
				,"gp1_eur_af":{"type":"float","index":"false"}
				,"gp1_af":{"type":"float","index":"true"}
				,"exac":{"type":"float","index":"true"}
				,"gmaf":{"type":"float","index":"false"}
				,"rd_freq":{"type":"float","index":"false"}
				,"gnomad_af":{"type":"float","index":"true"}
				,"gnomad_ac":{"type":"integer","index":"false"}
				,"gnomad_an":{"type":"integer","index":"false"}
				,"gnomad_af_popmax":{"type":"float","index":"false"}
				,"gnomad_ac_popmax":{"type":"integer","index":"false"}
				,"gnomad_an_popmax":{"type":"integer","index":"false"}
				,"gnomad_filter": {"type": "keyword"}
				,"clinvar_filter":{
					"type":"nested",
					"properties": {
						"clnsig":{"type":"keyword"}}}
				,"gene":{"type":"keyword"}
				,"transcript":{"type":"keyword"}
				,"protein_change":{"type":"keyword"}
				,"driver_statement":{"type":"keyword"}
				,"known_oncogenic_source":{"type":"keyword"}
				,"known_oncogenic_reference":{"type":"keyword"}
				,"onco_filter":{"type":"keyword"}
				,"consequence":{"type":"keyword"}
				,"effs":{
					"type":"nested",
						"properties":{
							"codon_change":{"type":"keyword","index":"false"}
							,"amino_acid_change":{"type":"keyword","index":"false"}
							,"amino_acid_length":{"type":"keyword","index":"false"}
							,"effect":{"type":"keyword"}
							,"effect_impact":{"type":"keyword"}
							,"exon_rank":{"type":"keyword","index":"false"}
							,"functional_class":{"type":"keyword","index":"false"}
							,"gene_coding":{"type":"keyword"}
							,"gene_name":{"type":"keyword"}
							,"transcript_biotype":{"type":"keyword"}
							,"transcript_id":{"type":"keyword"}}}
				,"samples_germline":{
					"type":"nested",
					"properties":{
						"dp":{"type":"float"}
						,"gq":{"type":"float"}
						,"ad":{"type":"keyword"}
						,"gt":{"type":"keyword"}
						,"sample":{"type":"keyword"}
						,"multi":{"type":"keyword","index":"false"}
						,"diploid":{"type":"keyword","index":"false"}}}
				,"samples_somatic":{
					"type":"nested",
					"properties":{
						"gt":{"type":"keyword"}
						,"dp_tumor":{"type":"float"}
						,"dp_control":{"type":"float"}
						,"ad_tumor":{"type":"keyword"}
						,"ad_control":{"type":"keyword"}
						,"sample":{"type":"keyword"}
						,"multi":{"type":"keyword","index":"false"}
						,"nprogs":{"type":"integer","index":"true"}
						,"progs":{"type":"keyword"}}}}}}}
	"""
	x = _create_index(host, port, index_name, data, user, pwd)

	if self is not None:
		if 'flags' not in vars(self):
			self.flags = {'index': (x == 200, x)}
		else:
			self.flags['index'] = (x == 200, x)
		return self
	else:
		return x

