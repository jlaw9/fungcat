# file to contain settings used by multiple scripts

import utils.file_utils as utils
import socket

# Ouputs Structure:
# ouputs - version1 - all - method 1 - results.txt
#       |          |      - method 2 - results.txt
#       |           - species1
#       |
#        - version2

ALLOWEDVERSIONS = [
    "2017_10-seq-sim",
    "2017_10-string",
    "2017_10-seq-sim-string",
    "2017_10-seq-sim-14",
    "2017_10-seq-sim-x5-string-14",
    "2017_10-seq-sim-x5-string",
    "2017_10-string-14",
    ]

SEQ_SIM_NETWORKS = {
    "2017_10-seq-sim": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    "2017_10-seq-sim-14": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    "2017_10-seq-sim-string": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    "2017_10-seq-sim-x5-string-14": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    "2017_10-seq-sim-x5-string": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    }

NETWORK_VERSION_INPUTS = {
    "2017_10-seq-sim": ['SEQ_SIM'],
    "2017_10-seq-sim-14": ['SEQ_SIM'],
    "2017_10-seq-sim-string": ['SEQ_SIM', 'STRING'],
    "2017_10-seq-sim-x5-string-14": ['SEQ_SIM', 'STRING'],
    "2017_10-seq-sim-x5-string": ['SEQ_SIM', 'STRING'],
    "2017_10-string": ['STRING'],
    "2017_10-string-14": ['STRING'],
    }


VERSION_SELECTED_STRAINS = {
    "2017_10-seq-sim": 'inputs/selected-strains.txt',
    "2017_10-seq-sim-string": 'inputs/selected-strains.txt',
    "2017_10-seq-sim-x5-string": 'inputs/selected-strains.txt',
    "2017_10-string": 'inputs/selected-strains.txt',
    "2017_10-seq-sim-14": 'inputs/selected-strains-14.txt',
    "2017_10-seq-sim-x5-string-14": 'inputs/selected-strains-14.txt',
    "2017_10-string-14": 'inputs/selected-strains-14.txt',
}
SELECTED_STRAINS = "inputs/selected-strains.txt"

ALGORITHM_OPTIONS = {
    "local-ova": "--ova local",
    "local-ovn": "--ovn local",
    #"fun-flow-ova": "--ova functional-flow",
    "fun-flow-ovn": "--ovn functional-flow",
    "sinksource-ova": "--ova sinksource",
    "sinksource-ovn": "--ovn sinksource",
    "genemania-ova": "--ova genemania",
}

#STRING_NETWORKS = {
#    "2017_10-string": ""
#    "2017_10-seq-sim-string": ""
#    }
# STRING networks are in the individual taxon dirs
#STRING_TAXON_UNIPROT = "%s/%%s/%%s-uniprot-links-v10.5-%%d.txt" % (STRING_TAXON_DIR)

NAME_TO_SHORTNAME = {
    "Neisseria gonorrhoeae FA 1090"                             : "Ng",
    "Peptoclostridium difficile / Clostridioides difficile 630" : "Cd",
    "Helicobacter pylori 85962"                                 : "Hp",
    "Klebsiella pneumoniae"                                     : "Kp",
    "Enterococcus faecium DO"                                   : "Ef",
    "Escherichia coli K-12"                                     : "Ec",
    "Haemophilus influenzae RD KW20"                            : "Hi",
    "Bordetella pertussis Tohama I"                             : "Bp",
    "Burkholderia cepacia"                                      : "Bc",
    "Clostridium botulinum A str. Hall"                         : "Cb",
    "Acinetobacter baumannii"                                   : "Ab",
    "Staphylococcus aureus"                                     : "Sa",
    "Vibrio cholerae O1 biovar El Tor str. N16961"              : "Bc",
    "Yersinia pestis"                                           : "Yp",
    "Streptococcus pyogenes"                                    : "Sp",
    "Pseudomonas aeruginosa"                                    : "Pa",
    "Salmonella typhimurium / Salmonella enterica"              : "Se",
    "Shigella dysenteriae serotype 1 (strain Sd197)"            : "Sd",
    "Mycobacterium tuberculosis"                                : "Mt",
    }
NAME_TO_TAX = {}
TAX_TO_NAME = {}

GOA_DIR = "/data/inputs/goa"
GOA_TAXON_DIR = "%s/taxon" % (GOA_DIR)

# input files
GO_FILE = "%s/2017-09-26-go.obo" % (GOA_DIR)
GOA_FILE = "%s/2017-09-26-goa_uniprot_all.gaf.gz" % (GOA_DIR)

# parsed input files
# for example: inputs/goa/taxon/22839/22839-goa.gaf
GOA_TAXON_FILE = "%s/%%s/%%s-goa.gaf" % (GOA_TAXON_DIR)
# file containing annotations propogated up the GO DAG hierarchy
# also formatted to be used as input for GAIN
FUN_FILE = "%s/%%s/%%s-goa-all-fun.txt" % (GOA_TAXON_DIR)
# file containing all annotations for the 19 species
GOA_ALL_FUN_FILE = "%s/all-taxon-goa.txt" % (GOA_TAXON_DIR)

# STRING directories and file templates
STRING_DIR = "/data/inputs/string"
STRING_TAXON_DIR = "%s/taxon" % (STRING_DIR)
STRING_FILE = "%s/protein.links.full.v10.5.txt.gz" % (STRING_DIR)
STRING_TO_UNIPROT = "%s/full_uniprot_2_string.04_2015.tsv" % (STRING_DIR)
# cutoff to be used for the STRING interactions
# Ranges from 150-1000
# 400 is considered a "Medium" cutoff
STRING_CUTOFF = 400

# Template for a species/taxon STRING file
# Last number is the cutoff used on interactions in this file
# for example: inputs/string/taxon/9606/9606.links.full.v10.5-400.txt
STRING_TAXON_FILE = "%s/%%s/%%s.links.full.v10.5-%%d.txt" % (STRING_TAXON_DIR)
# STRING FILE mapped to uniprot IDs
STRING_TAXON_UNIPROT = "%s/%%s/%%s-uniprot-links-v10.5-%%d.txt" % (STRING_TAXON_DIR)

# first column is uniprot ID, last is STRING ID
UNIPROT_TO_STRING = "inputs/protein-similarity/uniprot-species/uniprot-proteins-19-strains-plus-string.tab"
# first column is uniprot ID, second is organism ID
UNIPROT_TO_SPECIES = "inputs/protein-similarity/uniprot-species/2017-10-17-uniprot-prots-19-species.tab"

# code directory
BIORITHM_ubuntu16 = "/home/jeffl/src/c++/biorithm/trunk"
# use the ubuntu14 built version of GAIN on ubuntu 14 because the c++ packages are different
BIORITHM_ubuntu14 = "/home/jeffl/src/c++/biorithm-ubuntu14/trunk"

BIORITHM = BIORITHM_ubuntu16
ubuntu14_hosts = ['spode', 'cuthbert']
for host in ubuntu14_hosts:
    if host in socket.gethostname():
        BIORITHM = BIORITHM_ubuntu14

VERSION = ''
INPUTSPREFIX = ''
RESULTSPREFIX = ''
# processed network file for each version
NETWORK_template = "inputs/%s/%s-net.txt"


def set_version(version):
    global VERSION, INPUTSPREFIX, RESULTSPREFIX, NETWORK
    global NAME_TO_TAX, TAX_TO_NAME

    VERSION = version
    print "Using version '%s'" % (VERSION)

    INPUTSPREFIX = "inputs/%s" % VERSION
    RESULTSPREFIX = "outputs/%s" % VERSION
    NETWORK = NETWORK_template % (VERSION, VERSION)

    selected_strains = utils.readItemSet(VERSION_SELECTED_STRAINS[VERSION], 1)
    TAX_TO_NAME = utils.readDict(VERSION_SELECTED_STRAINS[VERSION], 1,2)
    NAME_TO_TAX = utils.readDict(VERSION_SELECTED_STRAINS[VERSION], 2,1)

    return INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains
