import json
import requests

def create_index(host,port,index_name,data,user,pwd):
    url = "http://" + host + ":" + port + "/" + index_name
    headers = {'Content-Type': 'application/json'}
    response = requests.put(url,data=data,headers=headers,auth=(user,pwd))
    print(response)


def create_index_snv(host,port,index_name,version,num_shards,num_replicas,user,pwd):
    data="""
          {"settings":{"index":{"number_of_shards":""" + num_shards + ""","number_of_replicas":""" + num_replicas + """, "refresh_interval":"1000ms"}}
            ,"mappings":{"""+"\"" + version + "\""+"""
            :{
            "properties":{
                "chrom":{"type":"integer","index":"true", "enabled":"false"}
                ,"pos":{"type":"integer","index":"true","enabled":"false"}
                ,"ref":{"type":"keyword","index":"false","enabled":"false"}
                ,"alt":{"type":"keyword","index":"false","enabled":"false"}
                ,"indel":{"type":"keyword","index":"true","enabled":"false"}
                ,"freqInt":{"type":"float","enabled":"false"}
                ,"rs":{"type":"keyword", "index":"true","enabled":"false"}
                ,"cadd_phred":{"type":"float","index":"true","enabled":"false"}
                ,"gerp_rs":{"type":"keyword","index":"false","enabled":"false"}
                ,"mt":{"type":"float","index":"false","enabled":"false"}
                ,"mutationtaster_pred":{"type":"keyword","enabled":"false"}
                ,"phylop46way_placental":{"type":"keyword","index":"false","enabled":"false"}
                ,"polyphen2_hvar_pred":{"type":"keyword","enabled":"false"}
                ,"polyphen2_hvar_score":{"type":"float","index":"false","enabled":"false"}
                ,"sift_pred":{"type":"keyword","enabled":"false"}
                ,"sift_score":{"type":"float","index":"false","enabled":"false"}
                ,"siphy_29way_pi":{"type":"keyword","index":"false","enabled":"false"}
                ,"UMD":{"type":"keyword"}
                ,"clinvar_clnsig":{"type":"keyword","index":"false","enabled":"false"}
                ,"clinvar_clnsigconf":{"type":"keyword","index":"false","enabled":"false"}
                ,"clinvar_id":{"type":"integer","index":"false","enabled":"false"}
                ,"gp1_afr_af":{"type":"float","index":"false","enabled":"false"}
                ,"gp1_asn_af":{"type":"float","index":"false","enabled":"false"}
                ,"gp1_eur_af":{"type":"float","index":"false","enabled":"false"}
                ,"gp1_af":{"type":"float","index":"true","enabled":"false"}
                ,"exac":{"type":"float","index":"true","enabled":"false"}
                ,"gmaf":{"type":"float","index":"false","enabled":"false"}
                ,"rd_freq":{"type":"float","index":"false","enabled":"false"}
                ,"gnomad_af":{"type":"float","index":"true","enabled":"false"}
                ,"gnomad_ac":{"type":"integer","index":"false","enabled":"false"}
                ,"gnomad_an":{"type":"integer","index":"false","enabled":"false"}
                ,"gnomad_af_popmax":{"type":"float","index":"false","enabled":"false"}
                ,"gnomad_ac_popmax":{"type":"integer","index":"false","enabled":"false"}
                ,"gnomad_an_popmax":{"type":"integer","index":"false","enabled":"false"}
                ,"gnomad_filter": {"type": "keyword","enabled":"false"}
                ,"clinvar_filter":{
                     "type":"nested",
                     "properties": {
                         "clnsig":{"type":"keyword"}}}
                ,"gene":{"type":"keyword"}
                ,"transcript":{"type":"keyword"}
                ,"protein_change":{"type":"keyword","enabled":"false"}
                ,"driver_statement":{"type":"keyword","enabled":"false"}
                ,"known_oncogenic_source":{"type":"keyword","enabled":"false"}
                ,"known_oncogenic_reference":{"type":"keyword","enabled":"false"}
                ,"onco_filter":{"type":"keyword"}
                ,"consequence":{"type":"keyword","enabled":"false"}
                ,"effs":{
                     "type":"nested",
                     "properties":{
                         "codon_change":{"type":"keyword","index":"false","enabled":"false"}
                         ,"amino_acid_change":{"type":"keyword","index":"false","enabled":"false"}
                         ,"amino_acid_length":{"type":"keyword","index":"false","enabled":"false"}
                         ,"effect":{"type":"keyword","enabled":"false"}
                         ,"effect_impact":{"type":"keyword","enabled":"false"}
                         ,"exon_rank":{"type":"keyword","index":"false","enabled":"false"}
                         ,"functional_class":{"type":"keyword","index":"false","enabled":"false"}
                         ,"gene_coding":{"type":"keyword","enabled":"false"}
                         ,"gene_name":{"type":"keyword","enabled":"false"}
                         ,"transcript_biotype":{"type":"keyword","enabled":"false"}
                         ,"transcript_id":{"type":"keyword"}}}
                ,"samples_germline":{
                     "type":"nested",
                     "properties":{
                         "dp":{"type":"float","enabled":"false"}
                         ,"gq":{"type":"float","enabled":"false"}
                         ,"ad":{"type":"keyword","enabled":"false"}
                         ,"gt":{"type":"keyword","enabled":"false"}
                         ,"sample":{"type":"keyword"}
                         ,"multi":{"type":"keyword","index":"false","enabled":"false"}
                         ,"diploid":{"type":"keyword","index":"false"}}}
                ,"samples_somatic":{
                     "type":"nested",
                     "properties":{
                         "gt":{"type":"keyword","enabled":"false"}
                         ,"dp_tumor":{"type":"float","enabled":"false"}
                         ,"dp_control":{"type":"float","enabled":"false"}
                         ,"ad_tumor":{"type":"keyword","enabled":"false"}
                         ,"ad_control":{"type":"keyword","enabled":"false"}
                         ,"sample":{"type":"keyword"}
                         ,"multi":{"type":"keyword","index":"false","enabled":"false"}
                         ,"nprogs":{"type":"integer","index":"true","enabled":"false"}
                         ,"progs":{"type":"keyword"}}}}}}}
    """
    create_index(host,port,index_name,data,user,pwd)


def create_index_cnv(host,port,index_name,version,num_shards,num_replicas,user,pwd):
    data="""
          {"settings":{"index":{"number_of_shards":""" + num_shards + ""","number_of_replicas":""" + num_replicas + ""","refresh_interval":"1000ms"}}
            ,"mappings":{"""+"\"" + version + "\""+"""
            :{
            "properties":{
                "chrom":{"type":"integer","index":"true"}
                ,"start":{"type":"integer","index":"true"}
                ,"end":{"type":"integer","index":"false"} 
                ,"type":{"type":"keyword","index":"false"}        
                ,"cnt":{"type":"integer","index":"true"}  
                ,"tool":{"type":"keyword","index":"true"}  
                ,"bf":{"type":"float","index":"true"}
                ,"DGV_goldstd_overlap":{"type":"keyword","index":"false"}
                ,"DGV_goldstd_coordinates":{"type":"keyword","index":"false"}
                ,"sample_id":{"type": "keyword","index":"true"}
                ,"omim_number":{"type":"integer", "index":"false"}
                ,"omim_phenotype":{"type":"keyword"}
                ,"reads_expected":{"type":"integer","index":"false"}
                ,"reads_observed":{"type":"integer","index":"false"}
                ,"reads_ratio":{"type":"float","index":"false"}
                ,"genes":{
                     "type":"nested",
                     "properties":{
                        "gene_name":{"type":"keyword"}
                      }}}}}}
    """
    create_index(host,port,index_name,data,user,pwd)
