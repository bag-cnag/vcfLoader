import json
import urllib2
from elasticsearch import Elasticsearch

def create_index(host,port,index_name,version):
    data="""
          {"settings":{"index":{"number_of_shards":8,"number_of_replicas":0,"refresh_interval":"1000ms"}}
            ,"mappings":{"""+"\"" + version + "\""+"""
            :{"_all":{"enabled":false},
            "properties":{"chrom":{"type":"integer","index":"true"},
            "pos":{"type":"integer","index":"true"}
            ,"ref":{"type":"keyword","index":"false"}
            ,"alt":{"type":"keyword","index":"false"}
            ,"indel":{"type":"keyword"}
            ,"freqInt":{"type":"keyword"}
            ,"effs":{"type":"nested","properties":{
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
            ,"transcript_id":{"type":"keyword"}
            }}
            ,"predictions":{"type":"nested",
            "properties":{"cadd_phred":{"type":"float","index":"true"}
            ,"gerp_rs":{"type":"keyword","index":"false"}
            ,"mt":{"type":"keyword","index":"false"}
            ,"mutationtaster_pred":{"type":"keyword"}
            ,"phylop46way_placental":{"type":"keyword","index":"false"}
            ,"polyphen2_hvar_pred":{"type":"keyword"}
            ,"polyphen2_hvar_score":{"type":"keyword","index":"false"}
            ,"sift_pred":{"type":"keyword"}
            ,"sift_score":{"type":"keyword","index":"false"}
            ,"siphy_29way_pi":{"type":"keyword","index":"false"}
            ,"UMD":{"type":"keyword"}
            ,"clinvar":{"type":"keyword","index":"false"}
            ,"clinvar_filter":{"type":"keyword"}
            ,"clnacc":{"type":"keyword","index":"false"},
            "rs":{"type":"keyword"}
            }},
            "populations":{"type":"nested",
            "properties":{"gp1_afr_af":{"type":"float","index":"false"}
            ,"gp1_asn_af":{"type":"float","index":"false"}
            ,"gp1_eur_af":{"type":"float","index":"false"}
            ,"gp1_af":{"type":"float","null_value":0.0}
            ,"esp6500_aa":{"type":"float","null_value":0.0}
            ,"esp6500_ea":{"type":"float","null_value":0.0}
            ,"exac":{"type":"float","null_value":0.0}
            ,"gmaf":{"type":"float","index":"false"}
            ,"rd_freq":{"type":"float","index":"false"}}}
            ,"samples":{"type":"nested",
            "properties":{"dp":{"type":"float"}
            ,"gq":{"type":"float"}
            ,"ad":{"type":"keyword"}
            ,"gt":{"type":"keyword"}
            ,"sample":{"type":"keyword"}
            ,"multi":{"type":"keyword","index":"false"},
            "diploid":{"type":"keyword","index":"false"}}}}}}}
          """
    url="http://"+host+":"+port+"/"+index_name
    es = Elasticsearch(hosts=[host])
    response = es.indices.create(index=index_name,ignore=400,body=data)
    print response
    
def delete_index(host,port,index_name,version):
    url="http://"+host+":"+port+"/"+index_name
    es = Elasticsearch(hosts=[host])
    es.indices.delete(index=index_name, ignore=[400, 404])
