def chromStrToInt(chrom):
    print("chromto int funciton")
    print(chrom)
    if (chrom =="MT") :
        return "23"
    elif (chrom=="X") :
        return "24"
    elif (chrom=="Y") :
        return "25"
    elif (chrom=="All"):
        print("ALL if condition")
        return ""
    else :
        return chrom
def oneFile(chrom):
    if (chrom=="All"):
        return ""
    else:
        return chrom
def update_version(url):
    splitted=url.split("/")
    old_version = splitted[len(splitted)-1]
    ver, rev = str(old_version).split('.')
    new_version=ver + '.' + str(int(rev)+1)
    return url.replace(old_version,new_version)


def version_bump( uri, increment = 'version' ):
    pos = { 'version': 0, 'revision': 1, 'iteration': 2 }[ increment.lower() ]

    splitted = uri.split('/')
    old_version = splitted[ len( splitted ) - 1 ]
    pack = str( old_version ).split( '.' )
    pack[ pos ] = str( int( pack[ pos ] ) + 1 )

    if increment.lower() == 'version':
        pack[ 1 ] = "0"

    if increment.lower() in ('version', 'revision'):
        pack[ 2 ] = "0"

    new_version = '.'.join( pack )
    
    return uri.replace( old_version, new_version )
    
def buildFileName(name, chrom):
    return name.replace("chromosome", chrom)

def buildDestinationVEP(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEP/"+fileName
    else:
        return destination+"/annotatedVEP_somatic/"+fileName

def buildDestinationdbNSFP(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEPdbnSFP/"+fileName
    else:
        return destination+"/annotatedVEPdbnSFP_somatic/"+fileName

def buildDestinationCADD(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEPdbnSFPCadd/"+fileName
    else:
        return destination+"/annotatedVEPdbnSFPCadd_somatic/"+fileName

def buildDestinationClinvar(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEPdbnSFPCaddClinvar/"+fileName
    else:
        return destination+"/annotatedVEPdbnSFPCaddClinvar_somatic/"+fileName

def buildDestinationGnomADEx(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad/"+fileName
    else:
        return destination+"/annotatedVEPdbnSFPCaddClinvarExGnomad_somatic/"+fileName

def buildDestinationExAC(destination, fileName, somatic = False):
    if not somatic:
        return destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadExAC/"+fileName
    else:
        return destination+"/annotatedVEPdbnSFPCaddClinvarExGnomadExAC_somatic/"+fileName

def buildDestinationTransform(destination, somatic = False):
    if not somatic:
        return destination+"/variants"
    else:
        return destination+"/variants_somatic"

def buildOriginToElastic(destination, chrom, somatic = False):
    if not somatic:
        return destination+"/variants/chrom="+chrom
    else:
        return destination+"/variants_somatic/chrom="+chrom

def buildOriginToElasticDenseMatrix(destination, nmtx, chrom, somatic = False):
    return destination+"/variants/chrom-" + chrom + "-mtx-" + str(nmtx)
