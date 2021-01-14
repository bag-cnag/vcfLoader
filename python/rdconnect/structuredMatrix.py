"""structuredMatrix

This module contains the functions used to create a sparse matrix and to append
experiments to an already existing sparse matrix.
"""


def append_to_sparse_matrix(self = None, config, log = VoidLog(), largeBatch = 500, smallBatch = 100):
    """ [...]


    process/moving_to
    combine/sparse_matrix_path

    'applications/datamanagement/ip'
    'applications/datamanagement/api_exp_status_list'

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
    largeBatch
    smallBatch
    """
    self, isConfig, isHl = check_class_and_config(None, config, hl, log)
    self.log.info('Entering step "append_to_sparse_matrix"')

    if not isConfig:
        self.log.error('No configuration was provided')
        raise NoConfigurationException('No configuration was provided')

    if not isHl:
        self.log.error('No pointer to HAIL module was provided')
        raise NoHailContextException('No pointer to HAIL module was provided')

    chrom = config['process/chrom']
    chrom = chrom_str_to_int(str(chrom))
    source_path = self.config['process/moving_to']
    sparse_path = self.config['combine/sparse_matrix_path']

    log.debug('> Argument "chrom" filled with "{}"'.format(chrom))
    log.debug('> Argument "source_path" filled with "{}"'.format(source_path))
    log.debug('> Argument "largeBatch" filled with "{}"'.format(largeBatch))
    log.debug('> Argument "smallBatch" filled with "{}"'.format(smallBatch))

    # Get experiments to load from DM

    url = config['applications/datamanagement/ip']
    if not url.startswith('http://') and not url.startswith('https://'):
        url = 'https://{0}'.format(url)

    url = config['applications/datamanagement/api_exp_status_list'].format(url)

    headers = { 
        'accept': 'application/json', 'Content-Type': 'application/json',
        'Authorization': 'Token {0}'.format(config['applications/datamanagement/token']),
        'Host': config['applications/datamanagement/host'] 
    }
    data = "{\"page\": 1, \"pageSize\": " + str(batch) + ", \"fields\": [\"RD_Connect_ID_Experiment\",\"mapping\",\"variantCalling\",\"genomicsdb\",\"hdfs\",\"es\",\"in_platform\"], \"sorted\":[{\"id\":\"RD_Connect_ID_Experiment\",\"desc\":false}], \"filtered\":[{\"id\":\"variantCalling\",\"value\":\"pass\"},{\"id\":\"rohs\",\"value\":\"pass\"},{\"id\":\"in_platform\",\"value\":\"waiting\"}]}"
    log.debug('> Querying DM using url "{0}"'.format(url))

    response = requests.post(url, data = data, headers = headers, verify = False)
    if response.status_code != 200:
        log.error('Query DM for experiment list resulted in a {} message'.format(str(response.status_code)))
        sys.exit(2)

    to_process = [ x['RD_Connect_ID_Experiment'] for x in json.loads(response.content)['items'] ]
    log.debug('> Obtained a total of "{}" samples to move'.format(len(to_process)))
    
    all_group = get.experiment_by_group(config, log, False)
    log.debug('> Obtained a total of "{}" samples for the group'.format(len(all_group)))
    
    to_process_group = [ x for x in all_group if x['RD_Connect_ID_Experiment'] in to_process ]

    print("hello! I got {} experiments to process from the total of {} I asked".format(len(to_process_group), largeBatch))


   print('-' * 25) 
   print(to_process[0])
   print('-' * 25)
   print(all_group[0])
   print('-' * 25)
   print(to_process_group[0])
   print('-' * 25)





def create_files_list( experiments, chrom, elastic_dataset ):
    """Creates a dictionary using RD-Connect Experiment ID as key and the its file as value."""
    prefix = 'hdfs://rdhdfs1:27000/test/rdconnect/gVCF'
    rst = {}
    for x in experiments:
        if x[ 'RD_Connect_ID_Experiment' ] not in rst.keys() and x[ 'elastic_dataset' ] == elastic_dataset:
            rst[ x[ 'RD_Connect_ID_Experiment' ] ] = prefix + '/' + x[ 'Owner' ] + "/" + x[ 'RD_Connect_ID_Experiment' ] + '/' + x[ 'RD_Connect_ID_Experiment' ] + '.' + chrom + '.g.vcf.bgz'
    return rst


def create_batches_sparse( list_of_ids, dict_of_paths, uri, smallSize = 100, largeSize = 1500 ):
    cnt = 0
    rst = []

    smallBatch = []
    largeBatch = []
    added = False
    bumpRev = False

    for idx, itm in enumerate( list_of_ids ):   
        if len( smallBatch ) >= smallSize:
            largeBatch.append( { 'uri': uri, 'batch': smallBatch } )
            cnt += smallSize
            smallBatch = []
            added = True

        if cnt >= largeSize:
            rst.append( { 'uri': uri, 'batches': largeBatch } )
            largeBatch = [ ]
            cnt = 0

        if added:
            if cnt + smallSize >= largeSize:
                uri = utils.version_bump( uri, 'revision' )
                bumpRev = True
            else:
                uri = utils.version_bump( uri, 'iteration' )
                bumpRev = False
            added = False
            
        smallBatch.append( { 'RD_Connect_ID_Experiment': itm[ 'RD_Connect_ID_Experiment' ],
            'Phenotips_ID': itm[ 'Phenotips_ID' ],
            'File': dict_of_paths[ itm[ 'RD_Connect_ID_Experiment' ] ]
        } )

    if len( smallBatch ) != 0:
        if not bumpRev:
            uri = utils.version_bump( uri, 'revision' )
        rst.append( { 'uri': uri, 'batches': [ { 'uri': uri, 'batch': smallBatch } ] } )
    return rst

def create_superbatches_sparse( list_of_uris ):
    rst = []
    first_uri = list_of_uris.pop( 0 )
    dst = utils.version_bump( first_uri, 'version' )
    for uri in list_of_uris:
        rst.append( { 'in_1': first_uri, 'in_2': uri, 'out': dst })
        first_uri = dst
        dst = utils.version_bump( dst, 'revision' )
    return rst


def combine_sparse_martix( uri_sm_1, uri_sm_2, destination_path ):
    print( '[combine_sparse_martix]: merging "{}" and "{}" and saving it to "{}"'.format( uri_sm_1, uri_sm_2, destination_path ) )
    sm_1 = hl.read_matrix_table( uri_sm_1 )
    sm_2 = hl.read_matrix_table( uri_sm_2 )
    comb = combine_gvcfs( [ sm_1 ] + [ sm_2 ] )
    comb.write( destination_path, overwrite = True )


def loadGvcf2( hl, experiments, destinationPath, gvcfStorePath, chrom, partitions ):
    #print("[loadGvcf] {} --> {}".format( str( len( experiments ) ), destinationPath ) )
    def transformFile( mt ):
        return transform_gvcf(mt.annotate_rows(
            info = mt.info.annotate( MQ_DP = hl.null( hl.tint32 ), VarDP = hl.null( hl.tint32 ), QUALapprox = hl.null( hl.tint32 ) )
        ))
    def importFiles( files ):
        x = hl.import_vcfs(
            files,
            partitions = interval[ 'interval' ], 
            reference_genome = interval[ 'reference_genome' ], 
            array_elements_required = interval[ 'array_elements_required' ]
        )
        return x

    interval = getIntervalByChrom( chrom, partitions )
    vcfs = [ transformFile( mt ) for mt in importFiles( [ x[ 'File' ] for x in experiments ] ) ]

    if gvcfStorePath == None:
        comb = combine_gvcfs( vcfs )
    else:
        gvcf_store = hl.read_matrix_table( gvcfStorePath )
        comb = combine_gvcfs( [ gvcf_store ] + vcfs )
    comb.write( destinationPath, overwrite = True )

def createSparseMatrix( group, url_project, host_project, token, prefix_hdfs, chrom, sz_small_batch, sz_large_batch, partitions_chromosome, gvcf_store_path, new_gvcf_store_path, gpap_id, gpap_token, is_playground ):
    """Iterates to create the sparse matrix."""
    #lgr = create_logger( 'createSparseMatrix', '' )
    if (new_gvcf_store_path is None or new_gvcf_store_path == '') and (gvcf_store_path is None or gvcf_store_path == ''):
        raise Exception('To properly run "createSparseMatrix" you have to provide the arguments "gvcf_store_path" or "new_gvcf_store_path".')


    # Get all the experiments that have to processed from data-management
    experiments_in_group = getExperimentByGroup( group, url_project, host_project, token, prefix_hdfs, chrom, sz_small_batch, is_playground )
    print('experiments_in_group', len( experiments_in_group ))
    print('\t', experiments_in_group[ : 2 ])
    experiment_status = getExperimentStatus( group, url_project, host_project, token, is_playground )
    print('experiment_status', len( experiment_status ))
    print('\t', experiment_status[ : 2 ])
    experiments_to_be_loaded = getExperimentsToProcess( experiment_status, experiments_in_group, check_hdfs = False )
    print('experiments_to_be_loaded', len( experiments_to_be_loaded ))
    print('\t', experiments_to_be_loaded[ : 2 ])

    # if is_playground:
    #     files_to_be_loaded = [ buildPathPlayground( prefix_hdfs, group, x[ 'RD_Connect_ID_Experiment' ], chrom ) for x in experiments_to_be_loaded ]
    # else:
    #     files_to_be_loaded = [ buildPath( prefix_hdfs, group, x[ 'RD_Connect_ID_Experiment' ], chrom ) for x in experiments_to_be_loaded ]
   
    # Having the multiple IDs (RD-Connect ID and PhenoTIPS/PhenoStore ID) we can create the path to the gVCF


    experiments_in_group = [ x for x in experiments_in_group if x[ 'elastic_dataset' ] ==  'rdcon_1488_670' ]
    files_to_be_loaded = create_files_list(experiments_in_group, str(chrom), "rdcon_1488_670")


    print('files_to_be_loaded', len( files_to_be_loaded.keys() ))
    print('\t', list( files_to_be_loaded.keys() )[ : 2 ])

    print('experiments_in_group (2)', len( experiments_in_group ))
    print('\t', experiments_in_group[ : 2 ])
    
    # The argument "new_gvcf_store_path" contains the path to the new sm that will be created from the blocks of 100 experiments and saved as 1k5
    # The argument "gvcf_store_path" will contain the last sm matrix that can be of any size and that will accumulate the old plus the new experiments

    list_of_batches = create_batches_sparse( experiments_in_group, files_to_be_loaded, new_gvcf_store_path, smallSize = sz_small_batch, largeSize = sz_large_batch )

    print('RUNNING STEP1 - CREATION OF CUMMULATIVE MATRICES OF {} EXPERIMENTS INCREMENTING {} EXPERIMENTS AT A TIME'.format( sz_large_batch, sz_small_batch ) )
    
    for idx, batch in enumerate( list_of_batches ):
        print(' > Processing large batch {}/{}'.format(idx, len( list_of_batches ) ) )
        # load each of the small batches of 100 experiments
        accum = None
        for idx, pack in enumerate( batch[ 'batches' ] ):
            print('     > Loading pack #{} of {} gVCF '.format( idx, len( pack[ 'batch' ] ) ) )
            for f in pack[ 'batch' ]:
                print(f)
            uri = '{}/chrom-{}'.format( pack[ 'uri' ], chrom )
            loadGvcf2( hl, pack[ 'batch' ], uri, accum, chrom, partitions_chromosome )
            accum = uri

    uris = [ b[ 'uri' ] for b in list_of_batches ]
    if not( gvcf_store_path is None or gvcf_store_path == '' ):
        uris = [ gvcf_store_path ] + uris

    print('RUNNING STEP2 - MERGING OF CUMMULATIVE MATRICES' )
    superbatches = create_superbatches_sparse( uris )
    for idx, pack in enumerate( superbatches ):
        combine_sparse_martix( '{}/chrom-{}'.format( pack[ 'in_1' ], chrom ), '{}/chrom-{}'.format( pack[ 'in_2' ], chrom ), '{}/chrom-{}'.format( pack[ 'out' ], chrom ) )
