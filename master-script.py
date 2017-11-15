#!/usr/bin/python

print("Importing libraries")

from optparse import OptionParser
import os
import sys
import gzip
from tqdm import tqdm
import src.utils.file_utils as utils
import src.fungcat_settings as f_settings
import src.fungcat_utils as f_utils
import glob
import networkx as nx
import shutil


# setup global variables
#GOA_DIR = f_settings.GOA_DIR


def main():
    global INPUTSPREFIX, RESULTSPREFIX, NETWORK
    global VERSION
    opts, args = parse_args()

    for version in opts.version:
        VERSION = version
        #print "Using version '%s'" % (VERSION)
        INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains = f_settings.set_version(VERSION)
        utils.checkDir(INPUTSPREFIX)

        #selected_strains = utils.readItemSet(opts.selected_strains, 1)
        print "%d selected strains: %s" % (len(selected_strains), ', '.join(selected_strains))

        # setup all of the GO annotations
        setup_inputs(selected_strains, forced=opts.forced)

        algorithms = opts.algorithm
        if opts.algorithm is None:
            # TODO setup options for algorithms
            #for algorithm in opts.algorithm:
            # Onevsnone local, onevsnone sink-source, functional flow onevsnone
            algorithms = [
                "local-ova",
                "local-ovn",
                "fun-flow-ovn",
                "sinksource-ova",
                "sinksource-ovn",
                "genemania-ova",
            ]
            #only_functions = "inputs/toxic-functions.txt"

        # now run GAIN on the selected_strains using a set of methods
        # run GAIN on each individual species if specified
        strains_to_run = None if opts.gain_per_species is False else selected_strains
        setup_and_run_gain(algorithms, opts.exp_name, selected_strains=strains_to_run, 
                           only_category=opts.only_category, only_functions=opts.only_functions,
                           ignore_evidence_codes=opts.ignore_evidence_code)

        if opts.viz_cv is True:
            plot_cv_results(algorithms, opts.exp_name)


def plot_cv_results(algorithms, exp_name):
    # combine all of the different algorithm's results, and then plot the precision recall 
    gain_results_dir = "%s/all" % (RESULTSPREFIX)
    out_dir = "%s/all/viz/%s" % (RESULTSPREFIX, exp_name)
    utils.checkDir(out_dir)
    gain_results_combined = "%s/db-%s-cv.txt" % (out_dir, exp_name)
    print "Combining all of the different algorithm's results to %s, and plotting the precision by recall" % (gain_results_combined)
    with open(gain_results_combined, 'w') as out:
        for algorithm in algorithms:
            gain_results_file = "%s/%s/pred/db-%s-cv.txt" % (gain_results_dir, algorithm, exp_name)
            print "\treading %s" % (gain_results_file)
            with open(gain_results_file, 'r') as f:
                for line in f:
                    out.write(line)

    # for now, use the eval-gain.pl perl script Murali wrote
    curr_dir = os.getcwd()
    os.chdir(out_dir)
    command = "perl %s/gain/misc/eval-gain.pl -cv-file db-%s-cv.txt" % (f_settings.BIORITHM_ubuntu16, exp_name)
    print command
    utils.runCommand(command)
    os.chdir(curr_dir)

    # now copy the plots to a more friendly viz directory
    out_dir2 = "outputs/viz/%s/%s" % (exp_name, VERSION)
    print "Copying plots to %s" % (out_dir2)
    utils.checkDir(out_dir2)
    for f in glob.glob("%s/*.pdf" % (out_dir)):
        shutil.copy(f, out_dir2)


def setup_network(selected_strains, forced=False):
    if os.path.isfile(NETWORK) and forced is False: 
        print "%s already exists. Use --forced to overwrite it" % (NETWORK)
        return
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
        print "Concatenate the STRING networks into a combined STRING network"
        print "\t%s" % (NETWORK)
        # concatenate the STRING networks into a combined STRING network
        # TODO figure out if/how we need to modify the weights of the STRING networks
        # so combining the seq sim network will be more effective
        for strain in tqdm(selected_strains):
            net_file = f_settings.STRING_TAXON_UNIPROT % (strain, strain, f_settings.STRING_CUTOFF)
            append_network_to_file(net_file, NETWORK)

        if 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
            print "Also adding the sequence similarity network"
            print "\t%s" % (f_settings.SEQ_SIM_NETWORKS[VERSION])
            print "\tlimiting to the %d selected strains" % (len(selected_strains))
            print "\tmultipyling the weights by 5 to get to a similar range as the STRING network"
            # concatenate the seq sim network to the end of the combind STRING network
            append_network_to_file(f_settings.SEQ_SIM_NETWORKS[VERSION], NETWORK, selected_strains=selected_strains, mod_weight=5)
    elif 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
        print "Setting up the sequence similarity network"
        # if this is only the sequence similarity network, then just setup a symbolic link
        if not os.path.isfile(NETWORK):
            if len(selected_strains) < 19:
                append_network_to_file(f_settings.SEQ_SIM_NETWORKS[VERSION], NETWORK, selected_strains=selected_strains)
            else:
                print "Adding a symbolic link from %s to %s" % (f_settings.SEQ_SIM_NETWORKS[VERSION],NETWORK)
                os.symlink(f_settings.SEQ_SIM_NETWORKS[VERSION],NETWORK)


def append_network_to_file(network_file, out_network_file, selected_strains=None, mod_weight=None):
    """
    *mod_weight*: multiply weight of all edges in the network_file by given amount
    """
    if selected_strains is not None:
        uniprot_to_species = f_utils.get_uniprot_to_species()
    #print "\t%s" % (f_settings.SEQ_SIM_NETWORKS[VERSION])
    #print "\tlimiting to the %d selected strains" % (len(selected_strains))
    # concatenate the seq sim network to the end of the combind STRING network
    lines_added = 0 
    with open(network_file, 'r') as f:
        with open(out_network_file, 'a') as out:
            for line in f:
                if selected_strains is None and mod_weight is None:
                    out.write(line)
                    lines_added += 1
                else:
                    a, b, w = line.rstrip().split('\t')[0:3]

                    if selected_strains is not None:
                        # Try limiting the interactions to those between species with STRING interactions
                        if uniprot_to_species[a] not in selected_strains or uniprot_to_species[b] not in selected_strains:
                            continue
                    if mod_weight is not None:
                        w = float(w)*mod_weight

                    out.write("%s\t%s\t%s\n" % (a, b, str(w)))
                    lines_added += 1

    print "\t%d interactions added" % (lines_added)
    return


def setup_and_run_gain(algorithms, exp_name, selected_strains=None, only_category=None, only_functions=None, ignore_evidence_codes=None):
    print "Running GAIN using these algorithms: %s" % (", ".join(algorithms))
    # first run GAIN using all species
    if selected_strains is None:
        out_dir = "%s/all" % (RESULTSPREFIX)
        # TODO instead of splitting the algorithms here and requiring GAIN to parse the network and annotations each time
        # I could write everything to the same file and then parse it after
        for algorithm in tqdm(algorithms):
            alg_out_dir = "%s/%s/pred" % (out_dir, algorithm)
            utils.checkDir(alg_out_dir)
            # TODO figure out a better way to name the experiment
            #exp_name = algorithm
            run_gain([algorithm], alg_out_dir, exp_name, NETWORK, 
                    f_settings.GO_FILE, f_settings.GOA_ALL_FUN_FILE,
                    only_category=only_category, only_functions=only_functions,
                    ignore_evidence_codes=ignore_evidence_codes)
    else:
        # Run GAIN on each strain individually using that strain's STRING network
        for strain in selected_strains:
            out_dir = "%s/%s" % (RESULTSPREFIX, strain)
            # TODO instead of splitting the algorithms here and requiring GAIN to parse the network and annotations each time
            # I could write everything to the same file and then parse it after
            for algorithm in tqdm(algorithms):
                alg_out_dir = "%s/%s/pred" % (out_dir, algorithm)
                utils.checkDir(alg_out_dir)
                #exp_name = algorithm
                run_gain([algorithm], alg_out_dir, exp_name, 
                        f_settings.STRING_TAXON_UNIPROT % (strain, strain, f_settings.STRING_CUTOFF), 
                        f_settings.GO_FILE, f_settings.FUN_FILE % (strain, strain),
                        only_category=only_category, only_functions=only_functions,
                        ignore_evidence_codes=ignore_evidence_codes)


def run_gain(algorithms, out_dir, exp_name, network, 
             go_file, goa_file, only_category=None, only_functions=None, ignore_evidence_codes=None, forced=False):
    out_file = "%s/db-%s-cv.txt" % (out_dir, exp_name)
    if os.path.isfile(out_file) and forced is False:
        print "\t%s already exists. (TODO) Use --foced to overwrite" % (out_file)
        return

    # now run GAIN
    command = "time %s/gain/gain " % (f_settings.BIORITHM) + \
              "  --go-file %s" % (go_file)+ \
              "  --interactions-file %s" % (network) + \
              "  --functions-file %s" % (goa_file) + \
              "  --output-directory %s" % (out_dir) + \
              "  --experiment-name %s" % (exp_name) + \
              "  --cross-validate --type weighted"

    if only_category is not None:
        command += "  --only-category %s" % (only_category)

    if only_functions is not None:
        command += "  --only-functions-file %s" % (only_functions)

    if ignore_evidence_codes is not None:
        for evidence_code in ignore_evidence_codes:
            command += "  --ignore-evidence-code %s" % (evidence_code)

    for algorithm in algorithms:
        command += "  %s" % (f_settings.ALGORITHM_OPTIONS[algorithm])

    utils.runCommand(command)
        

def setup_inputs(selected_strains, forced=False):
    """ Function to split GO annotations for the selected strains 
    And setup files for running GAIN
    """
    # setup all of the GO annotations
    split_annotations_per_taxon(selected_strains)
    transfer_annotations(selected_strains)
    combine_annotations(selected_strains, forced=forced)

    # TODO add mapping the STRING networks to UniProt
    # combine STRING networks for each species
    print f_settings.NETWORK_VERSION_INPUTS[VERSION]
    setup_network(selected_strains, forced=forced)

    # normalized_num_paths = (b-a) * (float(num_paths_per_edge[(u,v)] - min_num_paths) / float(max_num_paths - min_num_paths)) + a

    # setup the sequence similarity network (from Shiv)

    # now make sure the network is setup
    if not os.path.isfile(NETWORK):
        # TODO add symbolic link automatically
        # or combine sequence similarity and STRING networks
        print "Please setup network file %s" % (NETWORK)
        sys.exit()
        #if 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS:
            #setup_sequence_similarity_network()
    print "Using network: %s" % (NETWORK)
    return


#def setup_sequence_similarity_network():
#VERSION = "2017_10-seq-sim"
#net_file = f_settings.SEQ_SIM_NETWORKS[VERSION]
##G = nx.DiGraph()
##G = nx.read_weighted_edgelist(net_file, delimiter="\t", create_using=G)
#G = nx.DiGraph()
#lines = utils.readColumns(net_file, 1, 2, 3)
#G.add_weighted_edges_from(lines)
#print nx.info(G)
#G2 = nx.read_weighted_edgelist(net_file, delimiter="\t")
#print nx.info(G2)


def split_annotations_per_taxon(selected_strains):
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
                  "  --goa-annotations %s " % (f_settings.GOA_FILE) + \
                  "  --out-dir %s " % (f_settings.GOA_TAXON_DIR) + \
                  "  --selected-strains %s " % (strains_to_split_file) 
        utils.runCommand(command)


def transfer_annotations(selected_strains):
    # also propogate direct annotations up the DAG for each of the strains
    print "Propogating direct annotations up the DAG for each of the %d strains" % (len(selected_strains))
    #already_exist = True
    for strain in tqdm(selected_strains):
        out_file = f_settings.FUN_FILE % (strain, strain)
        if os.path.isfile(out_file):
            #print "\t%s already exists. Skipping" % (out_file)
            continue
        else:
            print "\twriting to %s" % (out_file)
            command = "%s/process-go/process-go " % (f_settings.BIORITHM) + \
                      "  --go-file %s " % (f_settings.GO_FILE) + \
                      "  --annotations-file %s " % (f_settings.GOA_TAXON_FILE % (strain, strain)) + \
                      "  --output-file %s " % (out_file) + \
                      "  --transitive-closure " + \
                      "  --direction-transitive-closure up "
            utils.runCommand(command)


def combine_annotations(selected_strains, out_file=f_settings.GOA_ALL_FUN_FILE, forced=False):
    # now combine all of the annotations into one file for all the species
    # TODO what if the user has added more strains? This file should include those 
    print "Combining annotations of all %d species to: %s" % (len(selected_strains), out_file)
    if os.path.isfile(out_file) and forced is False:
        print "\t%s already exists. Use --forced to overwrite" % (out_file)
        return

    with open(out_file, 'w') as out:
        write_header = True
        for strain in tqdm(selected_strains):
            with open(f_settings.FUN_FILE % (strain, strain), 'r') as f:
                header = f.readline()
                # write the header from the first file read
                if write_header:
                    out.write(header)
                    write_header = False
                for line in f:
                    out.write(line)


def parse_args():
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--version',type='string', action='append',
                      help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    parser.add_option('-i', '--goa-annotations', 
                      help="File containing goa annotations downloaded from UniProt in GAF format")
    parser.add_option('-o', '--out-dir', type='string', metavar='STR',
                      help="Directory to write the interactions of each species to.")
    #parser.add_option('-s', '--selected-strains', type='string', default=f_settings.SELECTED_STRAINS,
    #                  help="Uniprot reference proteome strains to perform analyses on. Default: %s" % (f_settings.SELECTED_STRAINS))
    parser.add_option('', '--gain-per-species', action="store_true", default=False,
                      help="Run GAIN on each individual species instead of all of them combined.")
    parser.add_option('', '--ignore-evidence-code', action="append",
                      help="GO evidence code to ignore. You may use this option multiple times.")
    parser.add_option('', '--only-category',  type='string',
                      help="Run GAIN using only the GO terms of the specified category. Valid options: 'f', 'p', 'c'")
    parser.add_option('', '--only-functions', type='string',
                      help="Run GAIN using only the functions in a specified file.")
    parser.add_option('', '--exp-name', type='string',
                      help="Experiment name to use when running GAIN.")
    parser.add_option('-a', '--algorithm', action="append", 
                      help="Algorithm(s) to use when running GAIN. If not specified, many algorithms will be used. Options: '%s'" % ("', '".join(f_settings.ALGORITHM_OPTIONS.keys())))
    parser.add_option('-v', '--viz-cv', action="store_true", default=False, 
                      help="Visulize the cross-validation results with precision/recall curves for each GO term")
    parser.add_option('', '--forced', action="store_true", default=False,
                      help="Force re-writing input and output files if they already exist")
    (opts, args) = parser.parse_args()

    if opts.exp_name is None:
        print "--exp-name required"
        sys.exit(1)

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print "ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS))
            sys.exit(1)

    return opts, args

if __name__ == "__main__":
    main()
