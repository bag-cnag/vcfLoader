import sys
import logging
import warnings
import hail as hl
import os,requests,json
from hail.experimental.vcf_combiner import *
from hail.experimental import full_outer_join_mt
from hail.experimental.vcf_combiner.vcf_combiner import combine_gvcfs
from hail.experimental.vcf_combiner.vcf_combiner import transform_gvcf
from rdconnect import utils
from rdconnect.annotations import truncateAt
from datetime import datetime
from subprocess import PIPE, Popen
from pyspark.sql import Row


def create_logger( name, path ):
    now = datetime.now()
    date_time = now.strftime("%y%m%d_%H%M%S")
    logger = logging.getLogger( name )
    logger.setLevel( logging.DEBUG )
    fh = logging.FileHandler( os.path.join( path, 'vcfLoader_{}_debug.log'.format( date_time ) ) )
    fh.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel( logging.DEBUG )
    formatter = logging.Formatter( '%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
    fh.setFormatter( formatter )
    ch.setFormatter( formatter )
    logger.addHandler( fh )
    logger.addHandler( ch )
    return logger


def resource(filename):
    return os.path.join(filename)


def getExperimentStatus( group, url_project, host_project, token, is_playground ):
    """Get the status information for all experiments allowed to be used by the token."""
    if not url_project.startswith( 'http://' ) and not url_project.startswith( 'https://' ):
        url_project = 'https://{0}'.format( url_project )
    if is_playground:
        url = "{0}/datamanagement_service/api/statusbyexperiment/?format=json&group={1}&user=dpiscia&owner=False".format( url_project, group )
    else:
        url = "{0}/datamanagement/api/statusbyexperiment/?format=json&group={1}&user=dpiscia&owner=False".format( url_project, group )
    headers = { 'Authorization': token, 'Host': host_project }
    print( 'getExperimentStatus: {0}'.format( url ) )
    resp = requests.get( url, headers = headers, verify = False )
    data = json.loads( resp.content )
    return data

def getExperimentByGroup( group, url_project, host_project, token, prefix_hdfs, chrom, max_items_batch, is_playground ):
    """Get all the experiments for a given group."""
    if not url_project.startswith( 'http://' ) and not url_project.startswith( 'https://' ):
        url_project = 'https://{0}'.format( url_project )
    if is_playground:
        url = "{0}/datamanagement_service/api/samplebygroup/?format=json&group={1}&user=dpiscia&owner=False&elastic_index=True&inplatform=False".format( url_project, group )
    else:
        url = "{0}/datamanagement/api/samplebygroup/?format=json&group={1}&user=dpiscia&owner=False&elastic_index=True&inplatform=False".format( url_project, group )
    headers = { 'Authorization': token, 'Host': host_project }
    print( 'getExperimentByGroup: {0}'.format( url ) )
    resp = requests.get (url, headers = headers, verify = False )
    data = json.loads( resp.content )
    return data


# def getExperimentByGroup( group, url_project, host_project, token, prefix_hdfs, chrom, max_items_batch ):
#     """Get all the experiments for a given group."""
#     if not url_project.startswith( 'http://' ) and not url_project.startswith( 'https://' ):
#         url_project = 'https://{0}'.format( url_project )
#     url = "{0}/datamanagement_service/api/samplebygroup/?format=json&group={1}&user=dpiscia&owner=False".format( url_project, group )
#     headers = { 'Authorization': token, 'Host': host_project }
#     print( 'getExperimentByGroup: {0}'.format( url ) )
#     resp = requests.get (url, headers = headers, verify = False )
#     data = json.loads( resp.content )
#     return data


def getExperimentsToProcess( experiment_status, experiment_available, check_hdfs = False ):
    """Given the experiments seen by the user as well as their status, returns the ones that are in HDFS and have to be processed."""
    experiment_status = [ x for x in experiment_status if x[ 'genomicsdb' ] == 'waiting' ]
    if check_hdfs:
        experiment_status = [ x for x in experiment_status if x[ 'hdfs' ] == 'pass' ]
    experiment_status_2 = [ x[ 'Experiment' ] for x in experiment_status ]
    experiment_available_2 = [ x[ 'RD_Connect_ID_Experiment' ] for x in experiment_available ]
    selected_experiments = [ x for x in experiment_available_2 if x in experiment_status_2 ]
    #print("Dropped experiments")
    #print([ x for x in experiment_available if x[ 'RD_Connect_ID_Experiment' ] not in selected_experiments ])
    #return [ x for x in experiment_available if x[ 'RD_Connect_ID_Experiment' ] in selected_experiments ]
    return experiment_available

# def create_files_list(experiments,chrom,elastic_dataset):
#     prefix="hdfs://rdhdfs1:27000/test/rdconnect/gVCF"
#     elastic_dataset="rdcon_1488_670"
#     return [ prefix+"/"+x['Owner']+"/"+x['RD_Connect_ID_Experiment']+'/'+x['RD_Connect_ID_Experiment']+'.'+chrom+'.g.vcf.bgz' for x in experiments if x[ 'elastic_dataset' ] == elastic_dataset ]


def create_files_list( experiments, chrom, elastic_dataset ):
    """Creates a dictionary using RD-Connect Experiment ID as key and the its file as value."""
    prefix = 'hdfs://rdhdfs1:27000/test/rdconnect/gVCF'
    rst = {}
    for x in experiments:
        if x[ 'RD_Connect_ID_Experiment' ] not in rst.keys() and x[ 'elastic_dataset' ] == elastic_dataset:
            rst[ x[ 'RD_Connect_ID_Experiment' ] ] = prefix + '/' + x[ 'Owner' ] + "/" + x[ 'RD_Connect_ID_Experiment' ] + '/' + x[ 'RD_Connect_ID_Experiment' ] + '.' + chrom + '.g.vcf.bgz'
    return rst

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


    experiments_in_group = [ x for x in experiments_in_group if x[ 'elastic_dataset' ] ==  'rdcon_2384_1' ]
    files_to_be_loaded = create_files_list(experiments_in_group, str(chrom), "rdcon_2384_1")


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
    

def create_batches_by_family( experiments, size = 1000 ):
    rst = []
    mtx = 0
    while len( experiments ) > 0:
        batch = []
        cnt = 0
        while cnt <= size and len( experiments ) > 0:
            fam = experiments[ 0 ][ 2 ]
            exp_fam = [ x for x in experiments if x[ 2 ] == fam ]
            for x in exp_fam:
                x.append('mtx' + str(mtx))
            batch += exp_fam
            cnt += len( exp_fam )
            experiments = [ x for x in experiments if x[ 2 ] != fam ]
        batch.sort(key = lambda x: x[ 0 ])
        rst += batch
        mtx += 1
    return rst


def create_family_groups(sc, sq, chrom, group, url_project, host_project, token, gpap_id,gpap_token,  prefix_hdfs, max_items_batch, sparse_matrix_path, dense_matrix_path, is_playground):
    lgr = create_logger('create_family_groups', '')
    chrom = "21"
    lgr.debug('OVERWRITING chrom to chrom-21')

    if sparse_matrix_path is None:
        raise 'No information on "sparse_matrix_path" was provided.'
    
    path_matrix = '{0}/chrom-{1}'.format( sparse_matrix_path, chrom )
    lgr.debug( 'READING from in {0}'.format( path_matrix ) )
    sparse_matrix = hl.read_matrix_table( path_matrix )
    
    experiments_in_matrix = [ x.get( 's' ) for x in sparse_matrix.col.collect() ]
    lgr.debug('Total of {0} experiments'.format( len( experiments_in_matrix ) ))

    # Get all the experiments that have to processed from data-management
    experiments_in_group = getExperimentByGroup( group, url_project, host_project, token, prefix_hdfs, chrom, max_items_batch, is_playground )
    print('experiments_in_group', len( experiments_in_group ))
    print('\t', experiments_in_group[ : 2 ])
    full_ids_in_matrix = [ x for x in experiments_in_group if x[ 'RD_Connect_ID_Experiment' ] in experiments_in_matrix ]
    print('full_ids_in_matrix', len( full_ids_in_matrix ))
    print('\t', full_ids_in_matrix[ : 2 ])
    experiments_and_families = getExperimentsByFamily( full_ids_in_matrix, url_project, gpap_id, gpap_token )
    print('experiments_and_families', len( experiments_and_families ))

    # Relocate experiments with no family
    none_detected = False
    x = len( list( set( [ x[ 2 ] for x in experiments_and_families ] ) ) )
    for ii in range( len( experiments_and_families ) ):
        if experiments_and_families[ ii ][ 2 ] == '---':
            none_detected = True
            experiments_and_families[ ii ][ 2 ] = experiments_and_families[ ii ][ 0 ]
    y = len( list( set( [ x[ 2 ] for x in experiments_and_families ] ) ) )
    if none_detected:
        warnings.warn( 'Provided experiment ids got no family assigned. RD-Connect ID used as family ID for those experiments. Original families were of {} while after update are of {}.'.format( x, y ) )
    experiments_and_families.sort(key=lambda x: x[ 0 ])

    batches = create_batches_by_family( experiments_and_families, 1000 )
    lgr.debug( 'Created {} batches'.format( len( batches ) ) )
    #for ii, bat in enumerate(batches):
    #    print('\tBatch {0}: {1} --> {2} - {3}'.format( ii, len( bat ), bat[0], bat[len(bat) - 1]))
    for sam in batches:
        print(sam[0], "\t", sam[1], "\t", sam[2], "\t", sam[3])

    log_path = '{0}/mapping'.format(dense_matrix_path)
    rdd = sc.parallelize(batches)
    experiments = rdd.map(lambda x: Row( RD_Connect_ID = x[ 0 ], PhenoTips = x[ 1 ], Family = x[ 2 ], DMatrix = x[ 3 ]))
    df = sq.createDataFrame(experiments)
    df.repartition(1).write.format('csv').mode('overwrite').save(log_path, header = 'true')


def load_table_log( sq, path ):
    sparlse_log = sq.read.format( 'csv' ).option( 'header', 'true' ).load( path )
    table = [ (str(row.RD_Connect_ID), str(row.PhenoTips), str(row.Family), int(str(row.DMatrix).replace("mtx", ""))) for row in sparlse_log.collect() ]
    mapping = []
    for mtx in list(set([ x[3] for x in table ])):
        y = [ x for x in table if x[3] == mtx ]
        #print("\t{} / {} : {} --> {}".format(mtx, len(y), y[0], y[len(y)- 1]))
        mapping.append(y)
    #print('table rows: {}'.format(len(table)))
    #print('mapping len: {}'.format(len(mapping)))
    mapping.sort(key=lambda x: x[0][3])
    return mapping



def createDenseMatrix( sc, sq, url_project, host_project, prefix_hdfs, dense_matrix_path, sparse_matrix_path, chrom, group, token, gpap_id, gpap_token, is_playground ):
    lgr = create_logger( 'createDenseMatrix', '' )

    mapping = load_table_log(sq, '{0}/mapping'.format(dense_matrix_path))

    if sparse_matrix_path is None:
        raise Exception('No information on "sparse_matrix_path" was provided.')
    
    path_matrix = '{0}/chrom-{1}'.format( sparse_matrix_path, chrom )
    lgr.debug( 'READING from in {0}'.format( path_matrix ) )
    sparse_matrix = hl.read_matrix_table( path_matrix )
    
    experiments_in_matrix = [ x.get( 's' ) for x in sparse_matrix.col.collect() ]
    lgr.debug('Total of {0} experiments'.format( len( experiments_in_matrix ) ))

    try:
        for idx, batch in enumerate( mapping ):
            lgr.debug( "Flatting and filtering dense matrix {0} (sz: {1}) --> {2} - {3}".format( idx, len( batch ), batch[0], batch[len(batch) - 1] ) )
            sam = hl.literal( [ x[ 0 ] for x in batch ], 'array<str>' )
            small_matrix = sparse_matrix.filter_cols( sam.contains( sparse_matrix['s'] ) )
            small_matrix = small_matrix.key_rows_by(small_matrix.locus, small_matrix.alleles)
            small_matrix = hl.experimental.sparse_split_multi( small_matrix, filter_changed_loci=True )
            small_matrix = hl.experimental.densify( small_matrix )
            small_matrix = small_matrix.filter_rows( hl.agg.any( small_matrix.GT.is_non_ref() ) )
            path = '{0}/chrom-{1}-mtx-{2}'.format( dense_matrix_path, chrom, idx )
            lgr.info( 'Writing dense matrix {} to disk ({})'.format( idx, path ) )
            small_matrix.write( path, overwrite = True )
            lgr.debug( "Ending writing dense matrix" )
    except Exception as ex:
        raise ex


def getExperimentsByFamily( pids, url_project, id_gpap, token_gpap, sort_output = True ):
    """Function to get the IDs from phenotips, from experiments, and from family."""
    print( "{0} ---> {1} / {2}".format( "getExperimentsByFamily", pids[ 0 ], pids[ len(pids) - 1 ] ) )
    #url = 'http://rdproto10:8082/phenotips/ExportMultiple'
    url = 'http://rdcompute3:8082/phenotips/ExportMultiple'
    data=[]
    headers = { 'Content-Type': 'application/json' }
    for i in range(0,(len(pids)//1000)+1) :
        body = { 'patients': [ { 'id': x[ 'Phenotips_ID' ] } for x in pids[(i*1000):((i+1)*1000)] ] }
        resp = requests.post( url, headers = headers, json = body, verify = False )
        data = data + resp.json()
    parsed = {}
    #import pdb;pdb.set_trace()
    for elm in data:
        pid = list( elm.keys() )[ 0 ]
        if type( elm[ pid ] ) == str:
            fam = '---'
        else: 
            fam = elm[ pid ][ 'family' ] if 'family' in elm[ pid ].keys() else '---'
        parsed[ pid ] = fam
    rst = [ [ pak[ 'RD_Connect_ID_Experiment' ], pak[ 'Phenotips_ID' ], parsed[ pak[ 'Phenotips_ID' ] ] ] for pak in pids ]
    if sort_output:
        return sorted( rst, key = lambda x: x[ 2 ] )
    else:
        return rst

def divideChunks( collection, size ): 
    # looping till length l 
    for ii in range( 0, len( collection ), size ):  
        yield collection[ ii:(ii + size) ] 


def get_samples(url_sample):
    splitted=url_sample.split("/")
    name_file=splitted[len(splitted)-1]
    return name_file.replace(".g.vcf.bgz","").split(".")


def getIntervals( chrom, max_pos, partitions ):
    quantity = max_pos // partitions
    intervals = []
    for item in range( 1, partitions + 1 ):
        start = (item - 1) * quantity
        end = item * quantity
        if item == len( range( 1, partitions + 1 ) ):
            end = max_pos
        if start == 0:
            start = 1
        intervals.append( hl.Interval( hl.Locus( str( chrom ), start ),hl.Locus( str( chrom ), end - 1 ), includes_end = True ) )
    return intervals


def getIntervalByChrom( chrom, partitions ):
    if chrom in ('23', 23):
        chrom = 'MT'
    if chrom in ('24', 24):
        chrom = 'X'
    if chrom in ('25', 25):
        chrom = 'Y'
    intervals = { # information from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37
        "25": { "interval": getIntervals( chrom,  59373566, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # Y
         "Y": { "interval": getIntervals( chrom,  59373566, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # Y
        "24": { "interval": getIntervals( chrom, 155270560, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # X
         "X": { "interval": getIntervals( chrom, 155270560, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # X
        "23": { "interval": getIntervals( chrom,     16570, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # MT
        "MT": { "interval": getIntervals( chrom,     16570, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }, # MT
        "22": { "interval": getIntervals( chrom,  51304566, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "21": { "interval": getIntervals( chrom,  48129895, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "20": { "interval": getIntervals( chrom,  63025520, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "19": { "interval": getIntervals( chrom,  59128983, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "18": { "interval": getIntervals( chrom,  78077248, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "17": { "interval": getIntervals( chrom,  81195210, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "16": { "interval": getIntervals( chrom,  90354753, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "15": { "interval": getIntervals( chrom, 102531392, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "14": { "interval": getIntervals( chrom, 107349540, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "13": { "interval": getIntervals( chrom, 115169878, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "12": { "interval": getIntervals( chrom, 133851895, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "11": { "interval": getIntervals( chrom, 135006516, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
        "10": { "interval": getIntervals( chrom, 135534747, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "9": { "interval": getIntervals( chrom, 141213431, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "8": { "interval": getIntervals( chrom, 146364022, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "7": { "interval": getIntervals( chrom, 159138663, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "6": { "interval": getIntervals( chrom, 171115067, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "5": { "interval": getIntervals( chrom, 180915260, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "4": { "interval": getIntervals( chrom, 191154276, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "3": { "interval": getIntervals( chrom, 198022430, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "2": { "interval": getIntervals( chrom, 243199373, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False },
         "1": { "interval": getIntervals( chrom, 249250621, partitions ), 'reference_genome': 'GRCh37', 'array_elements_required': False }
    }
    return intervals[ chrom ]


def loadGvcf( hl, files, chrom, destinationPath, gvcfStorePath, partitions, lgr ):
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
    lgr.debug( 'Got {} intervals for chrm {}'.format( len( interval ), chrom ) )

    lgr.debug( 'Importing {} files'.format( len( files ) ) )
    vcfs = [ transformFile( mt ) for mt in importFiles( files ) ]
    lgr.debug( 'Transformed files' )

    if gvcfStorePath == None:
        comb = combine_gvcfs( vcfs )
    else:
        gvcf_store = hl.read_matrix_table( gvcfStorePath )
        comb = combine_gvcfs( [ gvcf_store ] + vcfs )
    lgr.debug( 'Combined gVCF files' )
    lgr.debug( 'Saving sparse matrix to "{}"'.format( destinationPath ) )
    comb.write( destinationPath, overwrite = True )
    lgr.debug( 'Ended saving sparse matrix' )
    


#check if an experiment has been uploaded to hdfs
def buildPath( prefix, group, experiment, chrom):
    return '{0}/{1}/{2}/{3}'.format( prefix, group, experiment, utils.buildFileName( '{0}.chromosome.g.vcf.bgz'.format( experiment ), chrom ) )

def buildPathPlayground( prefix, group, experiment, chrom):
    if chrom in ('23', 23):
        chrom = 'MT'
    if chrom in ('24', 24):
        chrom = 'X'
    if chrom in ('25', 25):
        chrom = 'Y'
    return '{0}/{1}'.format( prefix, utils.buildFileName( '{0}.chromosome.g.vcf.bgz'.format( experiment ), chrom ) )

# def is_exp_uploaded(url_project,experiment,headers):
#     url=url_project+"/datamanagement_service/api/statusbyexperiment/?experiment="+experiment
#     resp=requests.get(url, headers=headers, verify=False)

#     if (resp.ok):
#         data=json.loads(resp.content)
#         if (data["hdfs"]=="pass"):
#             return True
        
#     else:
#         return False

# #get a list of file paths for a group
# def get_experiment_by_group(group,url_project,token,prefix_hdfs,chrom,max_items_batch):
#     headers = {'Authorization': token}
#     url=url_project+"/datamanagement_service/api/samplebygroup/?format=json&group="+group+"&user=dpiscia&owner=False"
#     resp=requests.get(url,headers=headers, verify=False)
#     data= json.loads(resp.content)
#     response=[]
#     counter=0
#     for exp in data:
#         if counter==max_items_batch:
#             break
#         if (is_exp_uploaded(url_project,exp["RD_Connect_ID_Experiment"],headers)):
#             counter=counter+1
#             response.append(build_path(prefix_hdfs,exp["Owner"],exp["RD_Connect_ID_Experiment"],chrom))
#     return response
