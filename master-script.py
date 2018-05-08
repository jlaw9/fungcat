#!/usr/bin/python

print("Importing libraries")

from optparse import OptionParser
import os
import sys
sys.path.append("src")
import gzip
from tqdm import tqdm
import src.utils.file_utils as utils
import src.fungcat_settings as f_settings
import src.fungcat_utils as f_utils
import glob
import networkx as nx
import shutil
import socket
import time


# setup global variables
#GOA_DIR = f_settings.GOA_DIR


def main():
    global INPUTSPREFIX, RESULTSPREFIX, NETWORK
    global VERSION
    opts, args = parse_args()
    start_time = time.time()

    for version in opts.version:
        VERSION = version
        #print("Using version '%s'" % (VERSION)
        INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains = f_settings.set_version(VERSION, epsilon=opts.convergence_cutoff)
        utils.checkDir(INPUTSPREFIX)

        #selected_strains = utils.readItemSet(opts.selected_strains, 1)
        print("%d selected strains: %s" % (len(selected_strains), ', '.join(selected_strains)))

        # setup all of the GO annotations
        setup_inputs(selected_strains, forced=opts.forced)

        if opts.exp_name is None:
            print("Finished setting up inputs\nSpecify --exp-name to run GAIN")
            sys.exit()

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
        # TODO should I have the only functions file be part of the exp_name
        #only_functions = "inputs/toxic-functions.txt"

        # now run GAIN on the selected_strains using a set of methods
        # run GAIN on each individual species if specified
        strains_to_run = None if opts.gain_per_species is False else selected_strains
        setup_and_run_gain(algorithms, opts.exp_name, selected_strains=strains_to_run,
                           only_category=opts.only_category, only_functions=opts.only_functions,
                           all_functions_file=opts.all_fun_file,
                           only_predictions=opts.only_predictions, 
                           only_cv=opts.only_cv, folds=opts.cross_validation_folds,
                           printonly=opts.printonly, forced=opts.forcealg)

    if opts.viz_cv is True:
        plot_cv_results(opts.version, algorithms, opts.exp_name, opts.only_functions)

    # total running time in hours 
    total_time = (time.time() - start_time)/(60.0*60)
    if opts.send_email is not None:
        wyatt = "wyatt" in socket.gethostname()
        subject = "FunGCAT Function Prediction"
        only_functions = utils.readItemSet(opts.only_functions, 1)
        message = "Finished running %d algorithm(s) on %d function(s) in %0.2fhr at `date`" %(len(algorithms), len(only_functions), total_time)
        # send an email saying the job finished
        email_command = f_utils.setupEmailCommand(opts.send_email, subject=subject, message=message, from_wyatt=wyatt)
        os.system(email_command)
        #utils.runCommand(email_command)


def plot_cv_results(versions, algorithms, exp_name, only_functions_file):
    """
    Call the script to plot the precision and recall for each specified GO term
    """
    import src.plot.plot_prec_rec as plot_prec_rec
    only_functions = utils.readItemSet(only_functions_file, 1)
    #only_functions = ["5488"]
    print("\nCalling main of src/plot/plot_prec_rec.py")
    plot_prec_rec.main(versions, exp_name, algorithms, only_functions, pdf=True, forced=True)


def plot_cv_results_perl(algorithms, exp_name):
    """
    Old function to call the perl precision/recall script
    """
    # combine all of the different algorithm's results, and then plot the precision recall
    gain_results_dir = "%s/all" % (RESULTSPREFIX)
    out_dir = "%s/all/viz/%s" % (RESULTSPREFIX, exp_name)
    utils.checkDir(out_dir)
    gain_results_combined = "%s/db-%s-cv.txt" % (out_dir, exp_name)
    print("Combining all of the different algorithm's results to %s, and plotting the precision by recall" % (gain_results_combined))
    with open(gain_results_combined, 'w') as out:
        for algorithm in algorithms:
            gain_results_file = "%s/%s/pred/db-%s-cv.txt" % (gain_results_dir, algorithm, exp_name)
            print("\treading %s" % (gain_results_file))
            with open(gain_results_file, 'r') as f:
                for line in f:
                    out.write(line)

    # for now, use the eval-gain.pl perl script Murali wrote
    curr_dir = os.getcwd()
    os.chdir(out_dir)
    command = "perl %s/gain/misc/eval-gain.pl -cv-file db-%s-cv.txt" % (f_settings.BIORITHM_ubuntu16, exp_name)
    print(command)
    utils.runCommand(command)
    os.chdir(curr_dir)

    # now copy the plots to a more friendly viz directory
    out_dir2 = "outputs/viz/%s/%s" % (exp_name, VERSION)
    print("Copying plots to %s" % (out_dir2))
    utils.checkDir(out_dir2)
    for f in glob.glob("%s/*.pdf" % (out_dir)):
        shutil.copy(f, out_dir2)


def setup_network(selected_strains, forced=False):
    # TODO ensure the sequence similarity network file contains only one copy of each edge
    # If not, take the lesser weight of the two directions
    if os.path.isfile(NETWORK) and forced is False:
        print("%s already exists. Use --forced to overwrite it" % (NETWORK))
        return
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
        print("Concatenate the STRING networks into a combined STRING network")
        print("\t%s" % (NETWORK))
        # concatenate the STRING networks into a combined STRING network
        # TODO figure out if/how we need to modify the weights of the STRING networks
        # so combining the seq sim network will be more effective
        for strain in tqdm(selected_strains):
            net_file = f_settings.STRING_TAXON_UNIPROT % (strain, strain, f_settings.STRING_CUTOFF)
            append_network_to_file(net_file, NETWORK)

        if 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
            print("Also adding the sequence similarity network")
            print("\t%s" % (f_settings.SEQ_SIM_NETWORKS[VERSION]))
            print("\tlimiting to the %d selected strains" % (len(selected_strains)))
            print("\tmultipyling the weights by 5 to get to a similar range as the STRING network")
            # concatenate the seq sim network to the end of the combind STRING network
            append_network_to_file(f_settings.SEQ_SIM_NETWORKS[VERSION], NETWORK, selected_strains=selected_strains, mod_weight=5)
    elif 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[VERSION]:
        print("Setting up the sequence similarity network")
        # if this is only the sequence similarity network, then just setup a symbolic link
        if not os.path.isfile(NETWORK):
            if len(selected_strains) < 19:
                append_network_to_file(f_settings.SEQ_SIM_NETWORKS[VERSION], NETWORK, selected_strains=selected_strains)
            else:
                print("Adding a symbolic link from %s to %s" % (f_settings.SEQ_SIM_NETWORKS[VERSION],NETWORK))
                os.symlink(f_settings.SEQ_SIM_NETWORKS[VERSION],NETWORK)


def append_network_to_file(network_file, out_network_file, selected_strains=None, mod_weight=None):
    """
    *mod_weight*: multiply weight of all edges in the network_file by given amount
    """
    if selected_strains is not None:
        uniprot_to_species = f_utils.get_uniprot_to_species()
    #print("\t%s" % (f_settings.SEQ_SIM_NETWORKS[VERSION])
    #print("\tlimiting to the %d selected strains" % (len(selected_strains))
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

    print("\t%d interactions added" % (lines_added))
    return


def setup_and_run_gain(algorithms, exp_name, selected_strains=None,
        only_category=None, only_functions=None, all_functions_file=f_settings.GOA_ALL_FUN_FILE, #ignore_evidence_codes=None,
        only_predictions=False, only_cv=False, folds=5, printonly=False, forced=False):
    print("Running GAIN using these algorithms: %s" % (", ".join(algorithms)))

    # UPDATE 2017-12-21: the ignore_evidence_codes functionality in GAIN does not appear to do anything.
    # TODO Ignore the evidence codes by creating a copy of the annotations file with those evidence codes removed
    # for now, I just created a file with IEA removed, but I should automate this 
#    if ignore_evidence_codes is None:
#        all_functions_file = f_settings.GOA_ALL_FUN_FILE
#    elif 'IEA' in ignore_evidence_codes and len(ignore_evidence_codes) == 1 and selected_strains is None:
#        all_functions_file = f_settings.GOA_ALL_FUN_FILE_NOIEA
#    elif 'EXP' in ignore_evidence_codes and len(ignore_evidence_codes) == 1 and selected_strains is None:
#        all_functions_file = f_settings.GOA_ALL_FUN_FILE_NOIEA.replace("-noiea","-comp")
#    elif 'COMP' in ignore_evidence_codes and len(ignore_evidence_codes) == 1 and selected_strains is None:
#        all_functions_file = f_settings.GOA_ALL_FUN_FILE_NOIEA.replace("-noiea","-exp")
#    elif len(ignore_evidence_codes) > 0:
#        print("\n--ignore-evidence-codes not yet implemented for codes besides IEA (and not implemented for --selected-strains). Quitting")
#        sys.exit()

    # first run GAIN using all species
    if selected_strains is None:
        out_dir = "%s/all" % (RESULTSPREFIX)
        # TODO instead of splitting the algorithms here and requiring GAIN to parse the network and annotations each time
        # I could write everything to the same file and then parse it after
        for algorithm in tqdm(algorithms):
            alg_out_dir = "%s/%s/%s" % (out_dir, algorithm, exp_name)
            utils.checkDir(alg_out_dir)
            # TODO figure out a better way to name the experiment
            #exp_name = algorithm
            run_gain([algorithm], alg_out_dir, exp_name, NETWORK,
                    f_settings.GO_FILE, all_functions_file,
                    only_category=only_category, only_functions=only_functions,
                    only_predictions=only_predictions, only_cv=only_cv,
                    folds=folds, printonly=printonly, forced=forced)
    else:
        # Run GAIN on each strain individually using that strain's STRING network
        for strain in selected_strains:
            out_dir = "%s/%s" % (RESULTSPREFIX, strain)
            # TODO instead of splitting the algorithms here and requiring GAIN to parse the network and annotations each time
            # I could write everything to the same file and then parse it after
            for algorithm in tqdm(algorithms):
                alg_out_dir = "%s/%s/%s" % (out_dir, algorithm, exp_name)
                utils.checkDir(alg_out_dir)
                #exp_name = algorithm
                run_gain([algorithm], alg_out_dir, exp_name,
                        f_settings.STRING_TAXON_UNIPROT % (strain, strain, f_settings.STRING_CUTOFF),
                        f_settings.GO_FILE, f_settings.FUN_FILE % (strain, strain),
                        only_category=only_category, only_functions=only_functions,
                        only_predictions=only_predictions, only_cv=only_cv,
                        folds=folds, printonly=printonly, forced=forced)


def run_gain(algorithms, out_dir, exp_name, network,
             go_file, goa_file, only_category=None, only_functions=None,
             only_predictions=False, only_cv=False, folds=5, printonly=False, forced=False):
    out_file = "%s/db-%s-cv.txt" % (out_dir, exp_name)
    if not forced and os.path.isfile(out_file) and os.path.getsize(out_file) > 50:
        tqdm.write("\t%s already exists. Use --forcealg to overwrite" % (out_file))
        return

    # now run GAIN
    command = "/usr/bin/time -v %s/gain/gain " % (f_settings.BIORITHM) + \
              "  --go-file %s" % (go_file)+ \
              "  --interactions-file %s" % (network) + \
              "  --functions-file %s" % (goa_file) + \
              "  --output-directory %s" % (out_dir) + \
              "  --experiment-name %s" % (exp_name) + \
              "  --type unweighted"
              # TODO UPDATE 2018-04-26: I thought I was using the weights, but I was not! GAIN expects them to be in the 4th column.
              # I need to either add an empty 3rd column also use the --weights-file option
              #"  --type weighted"

    if only_predictions is True:
        command += "  --only-predictions"
    if only_cv is True:
        command += "  --only-cv --cross-validate --cross-validate-fold %d" % (folds)
    else:
        # do both predictions and cross-validation
        command += "  --cross-validate --cross-validate-fold %d" % (folds)
    if only_category is not None:
        command += "  --only-category %s" % (only_category)

    if only_functions is not None:
        command += "  --only-functions-file %s" % (only_functions)

    # UPDATE 2017-12-21: the ignore_evidence_codes functionality in GAIN does not appear to do anything.
    #if ignore_evidence_codes is not None:
    #    for evidence_code in ignore_evidence_codes:
    #        command += "  --ignore-evidence-code %s" % (evidence_code)

    for algorithm in algorithms:
        command += "  %s" % (f_settings.ALGORITHM_OPTIONS[algorithm])
        # UPDATE 2018-02-15: set the weight of the artificial sink (negative) node for sinksource plus
        # The default of 1 is taking too long
        if algorithm == "sinksource-ovn":
            command += "  --ovn-sinksource-edge-weight %s" % (str(f_settings.VERSION_SS_LAMBDA[VERSION]))

    #command += "  &2> db-%s.log" % (exp_name)

    if printonly is True:
        print("Running: %s" % (command))
    else:
        utils.runCommand(command)


def setup_inputs(selected_strains, forced=False):
    """ Function to split GO annotations for the selected strains
    And setup files for running GAIN
    """
    # setup all of the GO annotations
    split_annotations_per_taxon(selected_strains)
    # TODO make an option for this somehow
    all_fun_file = f_settings.GOA_ALL_FUN_FILE_NOPARTOF
    if os.path.isfile(all_fun_file) and forced is False:
        print("\t%s already exists. Use --forced to overwrite" % (all_fun_file))
    else:
        transfer_annotations(selected_strains,
            out_file_template=f_settings.FUN_FILE_NOPARTOF, forced=forced)
        combine_annotations(selected_strains, out_file=all_fun_file)

    # TODO add mapping the STRING networks to UniProt
    # combine STRING networks for each species
    print(f_settings.NETWORK_VERSION_INPUTS[VERSION])
    setup_network(selected_strains, forced=forced)

    # normalized_num_paths = (b-a) * (float(num_paths_per_edge[(u,v)] - min_num_paths) / float(max_num_paths - min_num_paths)) + a

    # setup the sequence similarity network (from Shiv)

    # now make sure the network is setup
    if not os.path.isfile(NETWORK):
        # TODO add symbolic link automatically
        # or combine sequence similarity and STRING networks
        print("Please setup network file %s" % (NETWORK))
        sys.exit()
        #if 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS:
            #setup_sequence_similarity_network()
    print("Using network: %s" % (NETWORK))
    return


#def setup_sequence_similarity_network():
#VERSION = "2017_10-seq-sim"
#net_file = f_settings.SEQ_SIM_NETWORKS[VERSION]
##G = nx.DiGraph()
##G = nx.read_weighted_edgelist(net_file, delimiter="\t", create_using=G)
#G = nx.DiGraph()
#lines = utils.readColumns(net_file, 1, 2, 3)
#G.add_weighted_edges_from(lines)
#print(nx.info(G)
#G2 = nx.read_weighted_edgelist(net_file, delimiter="\t")
#print(nx.info(G2)


def split_annotations_per_taxon(selected_strains):
    # first make sure the goa annotations are parsed for the selected species
    strains_to_split = []
    for strain in selected_strains:
        strain_dir = "%s/%s" % (f_settings.GOA_TAXON_DIR, strain)
        if not os.path.isdir(strain_dir):
            #print("Splitting annotations for %s" % (strain_dir)
            strains_to_split.append(strain)

    if len(strains_to_split) > 0:
        print("Splitting GO annotations to each of the following %d strains: %s" % (len(strains_to_split), ', '.join(strains_to_split)))
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


def transfer_annotations(selected_strains, 
        go_file=f_settings.GO_FILE, gaf_taxon_file_template=f_settings.GOA_TAXON_FILE,
        out_file_template=f_settings.FUN_FILE, forced=False):
    """
    The templates should contain two '%s', one for the directory, and one for the file name
    """
    # also propogate direct annotations up the DAG for each of the strains
    print("Propogating direct annotations up the DAG for each of the %d strains" % (len(selected_strains)))
    #already_exist = True
    for strain in tqdm(selected_strains):
        out_file = out_file_template % (strain, strain)
        if forced is False and os.path.isfile(out_file):
            print("\t%s already exists. Skipping" % (out_file))
            continue
        else:
        #print("\twriting to %s" % (out_file))
            command = "%s/process-go/process-go " % (f_settings.BIORITHM) + \
                      "  --go-file %s " % (go_file) + \
                      "  --annotations-file %s " % (gaf_taxon_file_template % (strain, strain)) + \
                      "  --output-file %s " % (out_file) + \
                      "  --transitive-closure " + \
                      "  --direction-transitive-closure up "
            utils.runCommand(command)


def combine_annotations(selected_strains, out_file=f_settings.GOA_ALL_FUN_FILE):
    # now combine all of the annotations into one file for all the species
    # TODO what if the user has added more strains? This file should include those
    print("Combining annotations of all %d species to: %s" % (len(selected_strains), out_file))

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
    return


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
    parser.add_option('-e', '--convergence-cutoff', type='float', default=f_settings.DEFAULT_GAIN_EPSILON,
                      help="Convergence cutoff (epsilon) for SinkSource & GeneMANIA. Allowed options are: %s" % (', '.join(map(str, f_settings.ALLOWED_GAIN_EPSILON))))
    parser.add_option('', '--only-predictions', action="store_true", default=False,
                      help="Make predictions only")
    parser.add_option('', '--only-cv', action="store_true", default=False,
                      help="Perform cross-validation only")
    parser.add_option('-C', '--cross-validation-folds', type='int', default=5,
                      help="Perform cross validation using the specified # of folds. Default=5")
    #parser.add_option('', '--ignore-evidence-code', action="append",
    #                  help="GO evidence code to ignore. You may use this option multiple times.")
    parser.add_option('', '--all-fun-file', default=f_settings.GOA_ALL_FUN_FILE,
                      help="File containing all functions to use when running GAIN. This option replaced the --ignore-evidence-code option. Default=%s" % (f_settings.GOA_ALL_FUN_FILE))
    parser.add_option('', '--only-category', type='string',
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
                      help="Force re-writing input files if they already exist")
    parser.add_option('', '--forcealg', action="store_true", default=False,
                      help="Force re-running GAIN if the output files already exist")
    parser.add_option('', '--printonly', action="store_true", default=False,
                      help="Don't actually run GAIN. Just print the command to be run")
    parser.add_option('','--send-email', action='store', type='string',
                      help='Send an email (to specified recipient, from jeffreynlaw@gmail.com) after the script has finished')
    (opts, args) = parser.parse_args()

    if opts.version is None:
        print("--version required")
        sys.exit(1)

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)

    if opts.algorithm is not None:
        for alg in opts.algorithm:
            if alg not in f_settings.ALGORITHM_OPTIONS:
                print("ERROR: '%s' not an algorithm option. Algorithms are: %s." % (alg, ', '.join(f_settings.ALGORITHM_OPTIONS)))
                sys.exit(1)

    if opts.convergence_cutoff not in f_settings.ALLOWED_GAIN_EPSILON:
        print("ERROR: '%f' not a valid convergence cutoff value. Allowed settings: %s."
              % (opts.convergence_cutoff, ', '.join(map(str, f_settings.ALLOWED_GAIN_EPSILON))))
        sys.exit(1)

    return opts, args


if __name__ == "__main__":
    main()
