from rdconnect import utils

def annotationsVEP(hc,source, destination, vepPath):
    variants= hc.read(source)
    print("running vep")
    varAnnotated= variants.vep(vepPath)
    varAnnotated.write(destination,overwrite=True)
    #hc._jvm.core.annotations.vep(hc._jhc,"dd")
    #variants =


def dbnsfpTAble(hc,sourcePath,destinationPath):
    dbnsfpTable=hc.import_table(sourcePath+"/Users/dpiscia/spark/dbs/1.dbNSFP2.9.2a.2000lines.txt").annotate('variant= Variant(`#chr`,`pos(1-coor)`.toInt,`ref`,`alt`)').key_by('variant')
    dbnsfpTable.write(destinationPath,overwrite=True)


def annotatedbnsfp(variants, dbsfp_path):
    dbnsfp = hc.import_table("/Users/dpiscia/RD-repositories/data/output/dnsfp.kt")
    variants.annotate_variants_table(dbnfspVar,root='va.dbnsfp')