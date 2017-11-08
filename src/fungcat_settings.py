# file to contain settings used by multiple scripts

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
    ]

SEQ_SIM_NETWORKS = {
    "2017_10-seq-sim": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    "2017_10-seq-sim-string": "inputs/protein-similarity/2017-10-08-shiv-similarity-network-uniprot.txt",
    }

NETWORK_VERSION_INPUTS = {
        "2017_10-seq-sim": ['SEQ_SIM'],
        "2017_10-seq-sim-string": ['SEQ_SIM', 'STRING'],
        "2017_10-string": ['STRING'],
    }

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

TAX_TO_NAME = {
    "Neisseria gonorrhoeae FA 1090"                             : "242231",
    "Peptoclostridium difficile / Clostridioides difficile 630" : "272563",
    "Helicobacter pylori 85962"                                 :  "85962",
    "Klebsiella pneumoniae"                                     : "272620",
    "Enterococcus faecium DO"                                   : "333849",
    "Escherichia coli K-12"                                     :  "83333",
    "Haemophilus influenzae RD KW20"                            :  "71421",
    "Bordetella pertussis Tohama I"                             : "257313",
    "Burkholderia cepacia"                                      : "269482",
    "Clostridium botulinum A str. Hall"                         : "441771",
    "Acinetobacter baumannii"                                   : "509170",
    "Staphylococcus aureus"                                     :  "93061",
    "Vibrio cholerae O1 biovar El Tor str. N16961"              : "243277",
    "Yersinia pestis"                                           :    "632",
    "Streptococcus pyogenes"                                    : "301447",
    "Pseudomonas aeruginosa"                                    : "208964",
    "Salmonella typhimurium / Salmonella enterica"              :  "99287",
    "Shigella dysenteriae serotype 1 (strain Sd197)"            : "300267",
    "Mycobacterium tuberculosis"                                :  "83332",
    }

SELECTED_STRAINS = "inputs/selected-strains.txt"

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

# code directory
BIORITHM = "/home/jeffl/src/c++/biorithm/trunk"
BIORITHM_ubuntu14 = "/home/jeffl/src/c++/biorithm-ubuntu14/trunk"

VERSION = ''
INPUTSPREFIX = ''
RESULTSPREFIX = ''
# processed network file for each version
NETWORK = "inputs/%s/%s-net.txt"


def set_version(version):
    global VERSION, INPUTSPREFIX, RESULTSPREFIX, NETWORK

    VERSION = version
    print "Using version %s" % (VERSION)

    INPUTSPREFIX = "inputs/%s" % VERSION
    RESULTSPREFIX = "outputs/%s" % VERSION
    NETWORK = NETWORK % (VERSION, VERSION)

    return INPUTSPREFIX, RESULTSPREFIX, NETWORK
