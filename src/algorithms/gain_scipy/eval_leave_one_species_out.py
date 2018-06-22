
# script to see how well we can retrieve annotations of another species
# using the annotations of the other 18

from optparse import OptionParser,OptionGroup
from collections import defaultdict
import os
import sys
from tqdm import tqdm
import itertools
sys.path.append("src")
import utils.file_utils as utils
import fungcat_settings as f_settings
import run_algs
import alg_utils
import src.algorithms.setup_sparse_networks as setup
import igacat.go_term_prediction_examples.go_term_prediction_examples as go_examples
import algorithms.aptrank.run_birgrank as run_birgrank
from scipy import sparse
import numpy as np


def parse_args(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    # general parameters
    group = OptionGroup(parser, 'Main Options')
    group.add_option('','--version',type='string',
                     help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    group.add_option('-N','--net-file',type='string',
                     help="Network file to use. Can specify one per version. Default is the version's default network")
    group.add_option('-A', '--algorithm', action="append",
                     help="Algorithm for which to get predictions. Default is all of them. Options: '%s'" % ("', '".join(run_algs.ALGORITHMS)))
    group.add_option('', '--exp-name', type='string',
                     help="Experiment name to use when running GAIN.")
    group.add_option('', '--pos-neg-file', type='string', action='append',
                     help="File containing positive and negative examples for each GO term")
    group.add_option('', '--only-functions', type='string',
                     help="Run GAIN using only the functions in a specified file (should be in XX format i.e., without the leading GO:00).")
    group.add_option('-G', '--goterm', type='string', action="append",
                     help="Specify the GO terms to use (should be in GO:00XX format)")
    group.add_option('', '--non-pos-as-neg-eval', action='store_true', default=False,
                     help="By default, only the species' negatives are used when evaluating. " +
                     "This option will cause all non-positives to be treated as negative examples")
    parser.add_option_group(group)

    # parameters for running algorithms
    group = OptionGroup(parser, 'Algorithm options')
    group.add_option('', '--unweighted', action="store_true", default=False,
                     help="Option to ignore edge weights when running algorithms. Default=False (weighted)")
    group.add_option('-l', '--sinksourceplus-lambda', type=float, 
                     help="lambda parameter to specify the weight connecting the unknowns to the negative 'ground' node. Default=None")
    group.add_option('-k', '--k', type=int, action="append",
                     help="Top-k for Ripple, Default=200")
    group.add_option('-t', '--t', type=int, action="append",
                     help="t parameter for Ripple. Default=2")
    group.add_option('-s', '--s', type=int, action="append",
                     help="s parameter for Ripple. Default=200")
    group.add_option('-a', '--alpha', type=float, action="append",
                     help="Alpha insulation parameter. Default=0.8")
    group.add_option('', '--eps', type=float, action="append",
                     help="Stopping criteria for SinkSource")
    group.add_option('', '--max-iters', type=int, default=1000,
                     help="Maximum # of iterations for SinkSource. Default=1000")
    group.add_option('-e', '--epsUB', type=float, action="append",
                     help="Parameter to return the top-k if all other nodes have an UB - epsUB < the kth node's LB. Default=0")
    parser.add_option_group(group)

    # additional parameters
    group = OptionGroup(parser, 'Additional options')
    group.add_option('-b', '--obo-file', type='string', default=f_settings.GO_FILE,
                     help="GO OBO file which contains the GO DAG. Used if running AptRank/BirgRank. Default: %s" % (f_settings.GO_FILE))
    group.add_option('-T', '--taxon', type='string', action='append',
                      help="Specify the species taxonomy ID to use to be left out. Multiple may be specified. Otherwise, all species will be used")
    group.add_option('-W', '--num-pred-to-write', type='int', default=100,
                     help="Number of predictions to write to the file. If 0, none will be written. If -1, all will be written. Default=100")
    group.add_option('', '--only-pred', action="store_true", default=False,
                     help="Perform predictions only")
    group.add_option('', '--forcealg', action="store_true", default=False,
                     help="Force re-running algorithms if the output files already exist")
    group.add_option('', '--verbose', action="store_true", default=False,
                     help="Print additional info about running times and such")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args(args)

    if opts.exp_name is None or opts.pos_neg_file is None:
        print("--exp-name, --pos-neg-file, required")
        sys.exit(1)

    if opts.version is None:
        opts.version = ["2017_10-seq-sim"]

    if opts.algorithm is None:
        opts.algorithm = ['sinksource']

    if opts.version not in f_settings.ALLOWEDVERSIONS:
        print("ERROR: '%s' not an allowed version. Options are: %s." % (opts.version, ', '.join(f_settings.ALLOWEDVERSIONS)))
        sys.exit(1)

    #if opts.algorithm not in f_settings.ALGORITHM_OPTIONS:
    #    print "--algorithm %s is not a valid algorithm name" % (opts.algorithm)
    #    sys.exit(1)

    # TODO
    for alg in opts.algorithm:
        if alg not in run_algs.ALGORITHMS:
            print("ERROR: '%s' not a valid algorithm name. Algorithms are: '%s'." % (alg, ', '.join(run_algs.ALGORITHMS)))
            sys.exit(1)

    opts.k = opts.k if opts.k is not None else [200]
    opts.t = opts.t if opts.t is not None else [2]
    opts.s = opts.s if opts.s is not None else [200]
    opts.alpha = opts.alpha if opts.alpha is not None else [0.8]
    # default for deltaUBLB is None
    opts.eps = opts.eps if opts.eps is not None else [0.0001]
    opts.epsUB = opts.epsUB if opts.epsUB is not None else [0]

    return opts


def main(version, exp_name, W, prots, ann_matrix, goids,
         algorithms, opts, taxons=None, alpha=0.8, dag_matrix=None):
    # option to use negative examples when evaluating predictions
    #use_negatives_for_eval = True 
    # start with less GO terms
    #exp_name = "eval-species-%s-%d-%d-25iters" % (evidence_codes, cut1, cut2)
    if opts.non_pos_as_neg_eval is False:
        exp_name += "-use-neg" 
    else:
        exp_name += "-non-pos-neg" 

    # TODO make this more streamlined
    if 'birgrank' in algorithms:
        dag_matrix, pos_matrix, dag_goids = setup_h_ann_matrices(prots, opts.obo_file, opts.pos_neg_file, goterms=goterms)

    # set the version
    # selected species
    selected_species = utils.readDict(f_settings.VERSION_SELECTED_STRAINS[version], 1, 2)
    if taxons is None:
        taxons = selected_species.keys()

    # set the version for the UNIPROT_TO_SPECIES
    _, RESULTSPREFIX, _, _ = f_settings.set_version(version)
    print("Getting species of each prot from %s" % (f_settings.UNIPROT_TO_SPECIES))
    # for each of the 19 species, leave out their annotations 
    # and see how well we can retrieve them 
    uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1,2)
    # also build the reverse, but with the prot IDs
    node2idx = {n: i for i, n in enumerate(prots)}
    species_to_uniprot_idx = defaultdict(set)
    for p in uniprot_to_species:
        species_to_uniprot_idx[uniprot_to_species[p]].add(node2idx.get(p))

    print("Training/testing with %d species, %d goterms" % (len(taxons), len(goids)))

    for alg in algorithms:
        # TODO add other options to the output file 
        out_file = "%s/all/%s/%s/ground-truth-%sl%d-a%s.txt" % (
            RESULTSPREFIX, alg, exp_name,
            'unw-' if opts.unweighted else '', 
            0 if opts.sinksourceplus_lambda is None else opts.sinksourceplus_lambda,
            str(alpha).replace('.', '_'))

        if os.path.isfile(out_file) and opts.forcealg:
            if len(taxons) > 1: 
                print("Removing %s as results will be appended to it for each taxon" % (out_file))
                os.remove(out_file)
        elif os.path.isfile(out_file) and opts.only_pred is False:
            if len(taxons) > 1: 
                print("%s results file already exists. Use --forcealg to overwrite it.\nQuitting" % (out_file))
                sys.exit()
        elif opts.only_pred is False:
            print("Writing results to %s" % (out_file))

    # split the sets of positives/negatives into 19 sets with one species left out of each
    for s in tqdm(sorted(taxons)):
    #for s in selected_species:
        # leave this taxon out by removing its annotations
        # rather than a dictionary, build a matrix
        train_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
        test_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
        sp_goterms = []
        for i in range(len(goids)):
            pos, neg = alg_utils.get_goid_pos_neg(ann_matrix, i)
            pos = set(list(pos))
            neg = set(list(neg))
            # TODO I should limit these to the proteins in the network
            train_pos = pos - species_to_uniprot_idx[s]
            test_pos = pos & species_to_uniprot_idx[s]
            train_neg = neg - species_to_uniprot_idx[s]
            test_neg = neg & species_to_uniprot_idx[s]
            if len(train_pos) == 0 or len(test_pos) == 0 or \
               (len(train_neg) == 0 or len(test_neg) == 0):
                continue
            sp_goterms.append(i) 
            # build an array of the scores and set it in the goid sparse matrix of scores
            pos_neg_arr = np.zeros(len(prots))
            pos_neg_arr[list(train_pos)] = 1
            pos_neg_arr[list(train_neg)] = -1
            train_ann_mat[i] = pos_neg_arr
            pos_neg_arr = np.zeros(len(prots))
            pos_neg_arr[list(test_pos)] = 1
            pos_neg_arr[list(test_neg)] = -1
            test_ann_mat[i] = pos_neg_arr

        if 'birgrank' in algorithms:
            # TODO the matrix for BirgRank and the train_matrix do not have the same goids. 
            # I need to build a pos_mat with the train_matrix annotations
            train_pos_mat = sparse.lil_matrix(pos_matrix.shape)
            dag_goids2idx = {g: i for i, g in enumerate(dag_goids)}
            for i in range(len(goids)):
                dag_goid_idx = dag_goids2idx[goids[i]]
                train_pos_mat[dag_goid_idx] = train_ann_mat[i]
            # now set the negatives to 0 as birgrank doesn't use negatives
            train_pos_mat[train_pos_mat < 0] = 0

        #s_goterms = set(goid_pos.keys()) & set(test_goid_pos.keys())
        tqdm.write("\n" + "-"*30)
        #print("\n" + "-"*30)
        #print("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
        tqdm.write("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
            s, selected_species[s], len(sp_goterms), len(goids)))

        if len(sp_goterms) == 0:
            print("\tskipping")
            continue

        if opts.only_pred:
            print("Making predictions only and writing %d to a file" % (opts.num_pred_to_write))
            test_ann_mat = None
        elif opts.non_pos_as_neg_eval is True: 
            print("Evaluating using all non-ground-truth positives as false positives")
        else:
            print("Evaluating using only the ground-truth negatives predicted as positives as false positives")

        for alg in algorithms:
            # for now, use most of the defaults
            # change to running sinksource with 25 iterations
            # leave alpha at default 0.8
            alg_runner = run_algs.Alg_Runner(
                version, exp_name, W, prots, train_ann_mat, goids,
                algorithms=[alg], unweighted=opts.unweighted,
                ss_lambda=opts.sinksourceplus_lambda,
                eps_list=opts.eps, epsUB_list=opts.epsUB, max_iters=opts.max_iters,
                k_list=opts.k, t_list=opts.t, s_list=opts.s, a_list=[alpha],
                num_pred_to_write=opts.num_pred_to_write, verbose=opts.verbose, 
                forcealg=opts.forcealg)
            if alg == 'birgrank':
                goid_scores = alg_runner.run_aptrank_with_params(
                    train_pos_mat, dag_matrix, alg=alg, alpha=alpha) 
                curr_goids = dag_goids.copy() 
            else:
                curr_goids = goids.copy()
                goid_scores = alg_runner.main()
            # this will write an file containing the fmax for each goterm 
            # with the taxon name in the name of the file
            write_prec_rec = False 
            if len(taxons) == 1:
                print("Also writing prec/rec stats")
                write_prec_rec = True 

            # now evaluate 
            out_dir = "outputs/%s/all/%s/%s" % (version, alg, exp_name)
            utils.checkDir(out_dir)
            out_pref = "%s/ground-truth-%sl%d-" % (
                out_dir, 'unw-' if opts.unweighted else '',
                0 if alg_runner.ss_lambda is None else int(alg_runner.ss_lambda))
            alg_runner.evaluate_ground_truth(
                goid_scores, curr_goids, test_ann_mat, out_pref,
                non_pos_as_neg_eval=opts.non_pos_as_neg_eval,
                taxon=s, write_prec_rec=write_prec_rec)


def setup_h_ann_matrices(prots, obo_file, pos_neg_files, goterms=None):
    # parse the go_dags first as it also sets up the goid_to_category dictionary
    # TODO store the go dags as a file 
    go_dags = go_examples.parse_obo_file_and_build_dags(obo_file)

    # combine the matrices from bp and mf
    # would it make a difference running them together vs separately? I guess potentially it could
    # especially if I included part_of edges
    # TODO the hstack doesn't actually work yet
    dag_matrix = sparse.csr_matrix((0,0))
    ann_matrix = sparse.csr_matrix((0,0))
    goids = []
    # TODO build a matrix with the direct annotations (i.e., from the gaf file)
        # propagate the predictions(?)
    # for now, just use all of the propagated annotations
    # and then evaluate using the scores
    for pos_neg_file in pos_neg_files:
        if 'bp' in pos_neg_file:
            h = 'bp'
        elif 'mf' in pos_neg_file:
            h = 'mf'
        elif 'cc' in pos_neg_file:
            h = 'cc'
        curr_dag_matrix, curr_ann_matrix, curr_goids = run_birgrank.build_h_ann_matrices(
            prots, go_dags, pos_neg_files=[pos_neg_file], h=h, goterms=goterms)
        dag_matrix = sparse.hstack([dag_matrix,curr_dag_matrix])
        ann_matrix = sparse.hstack([ann_matrix,curr_ann_matrix])
        goids += curr_goids

    return dag_matrix, ann_matrix, goids


if __name__ == "__main__":
    opts = parse_args(sys.argv)
    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 

    #goid_pos, goid_neg = alg_utils.parse_pos_neg_files(opts.pos_neg_file, goterms=goterms) 
    # load the network matrix and protein IDs
    net_file = opts.net_file
    if net_file is None:
        _, _, net_file, _ = f_settings.set_version(opts.version) 
    W, prots = alg_utils.setup_sparse_network(net_file)

    # now build the annotation matrix
    ann_matrix, goids = setup.setup_sparse_annotations(opts.pos_neg_file, goterms, prots,
                            selected_species=None, taxon=None)

    for alpha in opts.alpha:
        main(opts.version, opts.exp_name, W, prots, ann_matrix, goids,
             opts.algorithm, opts, taxons=opts.taxon, alpha=alpha)
