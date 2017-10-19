#!/usr/bin/python

print("Importing libraries")

from optparse import OptionParser
import os
import sys
import gzip
from tqdm import tqdm
import src.utils.file_utils as utils
import src.fungcat_settings as f_settings


# setup global variables
#GOA_DIR = f_settings.GOA_DIR


def main():
    opts, args = parse_args()

    selected_strains = utils.readItemSet(opts.selected_strains, 1)

    setup_inputs(selected_strains)

    # STRING networks for each species

    # setup the sequence similarity network (from Shiv)

    # now run GAIN on the selected_strains using a set of methods


def setup_inputs(selected_strains):
    """ Function to split GO annotations for the selected strains 
    And setup files for running GAIN
    """
    # first make sure the goa annotations are parsed for the selected species
    strains_to_split = []
    for strain in selected_strains:
        strain_dir = "%s/%s" % (f_settings.GOA_TAXON_DIR, strain)
        if not os.path.isdir(strain_dir):
            #print "Splitting annotations for %s" % (strain_dir)
            strains_to_split.append(strain)

    if len(strains_to_split) > 0:
        print "Splitting GO annotations to each of the following %d strains: %s" % (len(strains_to_split), ', '.join(strains_to_split))
        utils.checkDir(f_settings.GOA_TAXON_DIR)
        strains_to_split_file = "%s/strains_to_split.txt" % (f_settings.GOA_TAXON_DIR)
        with open(strains_to_split_file, 'w') as out:
            out.write('\n'.join(strains_to_split))

        # split the annotations to a GAF file for the specified taxon/species 
        command = "python src/goa-split-to-species.py " + \
                " --goa-annotations %s " % (f_settings.GOA_FILE) + \
                " --out-dir %s " % (f_settings.GOA_TAXON_DIR) + \
                " --selected-strains %s " % (strains_to_split_file) 
        utils.runCommand(command)

    # also propogate direct annotations up the DAG for each of the strains
    print "Propogating direct annotations up the DAG for each of the %d strains" % (len(strains_to_split))
    #for strain in tqdm(strains_to_split):
    for strain in tqdm(selected_strains):
        command = "%s/process-go/process-go " % (f_settings.BIORITHM) + \
                " --go-file %s " % (f_settings.GO_FILE) + \
                " --annotations-file %s " % (f_settings.GOA_TAXON_FILE % (strain, strain)) + \
                " --output-file %s " % (f_settings.FUN_FILE % (strain, strain)) + \
                " --transitive-closure " + \
                " --direction-transitive-closure up "
        utils.runCommand(command)


def parse_args():
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--goa-annotations', 
                      help="File containing goa annotations downloaded from UniProt in GAF format")
    parser.add_option('-o', '--out-dir', type='string', metavar='STR',
                      help="Directory to write the interactions of each species to.")
    parser.add_option('-s', '--selected-strains', type='string', default='inputs/selected-strains.txt',
                      help="Uniprot reference proteome strains to perform analyses on")
    parser.add_option('', '--forced', action="store_true", default=False,
                      help="Force re-writing input and output files if they already exist")
    (opts, args) = parser.parse_args()

    return opts, args

if __name__ == "__main__":
    main()
