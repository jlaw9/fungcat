
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
import algorithms.setup_sparse_networks as setup
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
    #group.add_option('-N','--net-file',type='string',
    #                 help="Network file to use. Can specify one per version. Default is the version's default network")
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
    parser.add_option_group(group)

    # parameters for evaluating the results
    group = OptionGroup(parser, 'Evaluation options')
    group.add_option('', '--pos-neg-file-eval', type='string', action='append',
                     help="File containing positive and negative examples for each GO term used to evaluate for each species")
    group.add_option('', '--non-pos-as-neg-eval', action='store_true', default=False,
                     help="By default, only the species' negatives are used when evaluating. " +
                     "This option will treat all (non-positive) prots of the left-out species as negative examples")
    group.add_option('', '--keep-ann', action='store_true', default=False,
                     help="Don't leave out annotations when running the algorithms" +
                     "TODO allow for the option to run and evaluate everything together")
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
    # birgrank options
    group.add_option('', '--theta', type=float, default=.5,
                     help="BirgRank parameter: (1-theta) percent of Rtrain used in seeding vectors")
    group.add_option('', '--mu', type=float, default=.5,
                     help="BirgRank parameter: (1-mu) percent of random walkers diffuse from G via Rtrain to H")
    parser.add_option_group(group)

    # parameters for STRING networks
    group = OptionGroup(parser, 'STRING options (only used if STRING is part of the version)')
    #group.add_option('', '--only-combined', action="store_true", default=False,
    #        help="Use only the STRING combined network: \n\tcombined_score")
    #group.add_option('', '--only-core', action="store_true", default=False,
    #        help="Use only the 6 core networks: \n\t%s" % (', '.join(setup.CORE_STRING_NETWORKS)))
    #group.add_option('', '--non-transferred', action="store_true", default=False,
    #        help="Use all non-transferred networks: \n\t%s" % (', '.join(NON_TRANSFERRED_STRING_NETWORKS)))
    #group.add_option('', '--all-string', action="store_true", default=False,
    #        help="Use all individual 13 STRING networks: \n\t%s" % (', '.join(STRING_NETWORKS)))
    group.add_option('-S', '--string-networks', type='string', default=', '.join(setup.CORE_STRING_NETWORKS),
            help="Comma-separated list of string networks to use. " +
                 "If specified, other STRING options will be ignored." +
                 "Default: %s" % (', '.join(setup.CORE_STRING_NETWORKS)))
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
         algorithms, opts, taxons=None, alpha=0.8, dag_matrix=None,
         eval_ann_matrix=None):
    """
    *W*: If using STRING networks, then W should be a tuple containing
        a list of all STRING sparse matrices to use when weighting (SWSN)
        a list of the names of the STRING networks
    """
    # option to use negative examples when evaluating predictions
    #use_negatives_for_eval = True 
    # start with less GO terms
    #exp_name = "eval-species-%s-%d-%d-25iters" % (evidence_codes, cut1, cut2)
    if opts.non_pos_as_neg_eval is False:
        exp_name += "-use-neg" 
    else:
        exp_name += "-non-pos-neg" 
    if opts.keep_ann:
        exp_name += "-keep-ann" 

    # TODO make this more streamlined
    if 'birgrank' in algorithms:
        dag_matrix, pos_matrix, dag_goids = setup_h_ann_matrices(
                prots, opts.obo_file, opts.pos_neg_file, goterms=goids)
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
        sparse_networks, net_names = W

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
        tqdm.write("\n" + "-"*30)
        #print("\n" + "-"*30)
        #print("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
        tqdm.write("Taxon: %s - %s" % (
            s, selected_species[s]))
        # leave this taxon out by removing its annotations
        # rather than a dictionary, build a matrix
        # convert the annotation matrix to a lil matrix
        ann_matrix = ann_matrix.tolil()
        if eval_ann_matrix is not None:
            eval_ann_matrix = eval_ann_matrix.tolil()
        train_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
        test_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
        sp_goterms = []
        for i in range(len(goids)):
            pos, neg = alg_utils.get_goid_pos_neg(ann_matrix, i)
            ann_pos = set(list(pos))
            ann_neg = set(list(neg))
            # first setup the training annotations (those used as positives/negatives for the algorithm)
            if opts.keep_ann:
                train_pos = ann_pos 
                train_neg = ann_neg 
            else:
                train_pos = ann_pos - species_to_uniprot_idx[s]
                train_neg = ann_neg - species_to_uniprot_idx[s]
            eval_pos = ann_pos.copy()
            eval_neg = ann_neg.copy()
            # setup the testing annotations (those used when evaluating the performance
            if eval_ann_matrix is not None:
                eval_pos, eval_neg = alg_utils.get_goid_pos_neg(eval_ann_matrix, i)
                eval_pos = set(list(eval_pos))
                eval_neg = set(list(eval_neg))
            # TODO I should limit these to the proteins in the network
            test_pos = eval_pos & species_to_uniprot_idx[s]
            # UPDATE 2018-06-27: Only evaluate the species prots as negatives, not all prots
            if opts.non_pos_as_neg_eval:
                test_neg = species_to_uniprot_idx[s] - eval_pos
                test_neg.discard(None)
                #test_neg = np.asarray(sorted(test_neg)).astype(int)
            else:
                test_neg = eval_neg & species_to_uniprot_idx[s]
            # UPDATE 2018-06-30: Remove test positives/negatives that are part of the training positives/negatives
            # don't remove test positives if its a training negative because not all algorithms use negatives
            test_pos -= train_pos 
            test_neg -= train_pos | train_neg 
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

        tqdm.write("\t%d/%d goterms with > 0 annotations" % (len(sp_goterms), len(goids)))

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
            # not needed for lil matrix
            #train_pos_mat.eliminate_zeros() 

        if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
            # use the simultaneous weighting method to weight the networks
            out_file = "inputs/%s/%s/%d-nets-combined-SWSN-leave-out-%s.npz" % (
                version, exp_name, len(sparse_networks), s)
            if os.path.isfile(out_file):
                print("Loading SWSN weighted network from %s" % (out_file))
                W = sparse.load_npz(out_file)
            else:
                # remove rows with 0 annotations
                empty_rows = []
                for i in range(len(goids)):
                    pos, neg = alg_utils.get_goid_pos_neg(train_ann_mat, i)
                    # the combineWeightsSWSN method doesn't seem to
                    # work if there's only 1 positive
                    if len(pos) <= 1 or len(neg) <= 1:
                        empty_rows.append(i)
                # don't modify the original to keep the rows matching the GO ids
                curr_train_ann_mat = alg_utils.delete_rows_csr(train_ann_mat.tocsr(), empty_rows)
                utils.checkDir(os.path.dirname(out_file))
                W = setup.weight_SWSN(curr_train_ann_mat, sparse_networks,
                        net_names=net_names, out_file=out_file, nodes=prots)

        if len(sp_goterms) == 0:
            print("\tskipping")
            continue

        if opts.keep_ann:
            print("Keep all annotations when making predictions")
        if opts.only_pred:
            print("Making predictions only and writing %d to a file" % (opts.num_pred_to_write))
            test_ann_mat = None
        elif opts.non_pos_as_neg_eval is True: 
            print("Evaluating using all non-ground-truth positives for the taxon as false positives")
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
                forcealg=opts.forcealg, progress_bar=False)
            if alg == 'birgrank':
                # the W matrix is already normalized, so I can run
                # birgrank/aptrank from here
                goid_scores = alg_runner.run_aptrank_with_params(
                    train_pos_mat, dag_matrix, alg=alg, alpha=alpha,
                    theta=opts.theta, mu=opts.mu) 
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
            out_dir = "outputs/%s/all/%s/%s/goids/" % (version, alg, exp_name)
            utils.checkDir(out_dir)
            out_pref = "%s/ground-truth-%sl%d-" % (
                out_dir, 'unw-' if opts.unweighted else '',
                0 if alg_runner.ss_lambda is None else int(alg_runner.ss_lambda))
            if alg == 'birgrank':
                out_pref += 'a%s-t%s-m%s' % (
                    str(alpha).replace('.','_'), str(opts.theta).replace('.','_'),
                    str(opts.mu).replace('.','_'))
            alg_runner.evaluate_ground_truth(
                goid_scores, curr_goids, test_ann_mat, out_pref,
                #non_pos_as_neg_eval=opts.non_pos_as_neg_eval,
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


def run():
    opts = parse_args(sys.argv)
    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 

    if 'birgrank' in opts.algorithm and len(opts.pos_neg_file) > 1:
        print("Birgrank not yet implemented with multiple hierarchies. Use only bp or mf")
        sys.exit()

    #goid_pos, goid_neg = alg_utils.parse_pos_neg_files(opts.pos_neg_file, goterms=goterms) 
    # load the network matrix and protein IDs
    #net_file = opts.net_file
    #if net_file is None:
    INPUTSPREFIX, _, net_file, selected_strains = f_settings.set_version(opts.version) 
    # TODO this should be better organized so that any STRING networks
    # can be used
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
        string_nets = setup.CORE_STRING_NETWORKS
        out_pref_net = "%s/sparse-nets/" % (INPUTSPREFIX)
        utils.checkDir(out_pref_net)
        # build the file containing the sparse networks
        sparse_networks, network_names, prots = setup.create_sparse_net_file(
            opts.version, out_pref_net, selected_strains=selected_strains,
            string_nets=string_nets, string_cutoff=f_settings.STRING_CUTOFF,
            forcenet=False)
        # TODO organize this better
        W = (sparse_networks, network_names)
    else:
        W, prots = alg_utils.setup_sparse_network(net_file)

    # now build the annotation matrix
    ann_matrix, goids = setup.setup_sparse_annotations(opts.pos_neg_file, goterms, prots,
                            selected_species=None, taxon=None)
    #print(goids)
    test_ann_matrix = None 
    if opts.pos_neg_file_eval is not None:
        test_ann_matrix, test_goids = setup.setup_sparse_annotations(opts.pos_neg_file_eval, goterms, prots,
                                selected_species=None, taxon=None)
        # Make sure the goids are the same before using it
        num_ann_with_test = len(set(goids) & set(test_goids))
        if num_ann_with_test != 0:
            #print(test_goids)
            print("WARNING: %d / %d goids in the annotation matrix are in the evaluation matrix" % (num_ann_with_test, len(goids)))
            ann_goids_without_eval = set(goids) - set(test_goids)

            print("Removing %d GO terms from the annotation matrix" % (len(ann_goids_without_eval)))
            goids_to_remove = []
            new_goids = []
            for i, g in enumerate(goids):
                if g in ann_goids_without_eval:
                    goids_to_remove.append(i)
                    continue
                new_goids.append(g)
            ann_matrix = alg_utils.delete_rows_csr(ann_matrix.tocsr(), goids_to_remove)
            goids = new_goids

        test_goids2idx = {g: i for i, g in enumerate(test_goids)}
        # slim down the evaluation matrix to have the same GO terms as the ann matrix
        matching_test_matrix = sparse.lil_matrix(ann_matrix.shape)
        for i, g in enumerate(goids):
            matching_test_matrix[i] = test_ann_matrix[test_goids2idx[g]]
        test_ann_matrix = matching_test_matrix

    # TODO alpha doesn't change for all algorithms, just SS
    for alpha in opts.alpha:
        main(opts.version, opts.exp_name, W, prots, ann_matrix, goids,
             opts.algorithm, opts, taxons=opts.taxon, alpha=alpha,
             eval_ann_matrix=test_ann_matrix)


if __name__ == "__main__":
    run()
