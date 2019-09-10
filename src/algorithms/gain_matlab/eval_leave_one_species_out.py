
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
#import run_algs
import algorithms.gain_matlab.run_algs as run_algs
import algorithms.gain_matlab.alg_utils as alg_utils
import algorithms.setup_sparse_networks as setup
import algorithms.aptrank.run_birgrank as run_birgrank
from scipy import sparse
import numpy as np
import pdb


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
    group.add_option('', '--pos-neg-file', type='string',
                     help="File containing positive and negative examples for each GO term")
    group.add_option('', '--only-functions', type='string',
                     help="Run GAIN using only the functions in a specified file (should be in XX format i.e., without the leading GO:00).")
    group.add_option('-G', '--goterm', type='string', action="append",
                     help="Specify the GO terms to use (should be in GO:00XX format)")
    parser.add_option_group(group)

    # parameters for evaluating the results
    group = OptionGroup(parser, 'Evaluation options')
    group.add_option('', '--pos-neg-file-eval', type='string',
                     help="File containing positive and negative examples for each GO term used to evaluate for each species")
    group.add_option('', '--non-pos-as-neg-eval', action='store_true', default=False,
                     help="By default, only the species' negatives are used when evaluating. " +
                     "This option will treat all (non-positive) prots of the left-out species as negative examples")
    group.add_option('', '--oracle', action='store_true', default=False,
                     help="option to remove train negatives that are actually test positives")
    group.add_option('', '--keep-ann', action='store_true', default=False,
                     help="Don't leave out annotations when running the algorithms" +
                     "TODO allow for the option to run and evaluate everything together")
    group.add_option('', '--num-test-cutoff', type='int', default=1,
                     help="Minimum number of annotations for each GO term in the left-out species to test. Default: 1")
    group.add_option('', '--eval-goterms-with-left-out-only', action='store_true', default=False,
                     help="if --pos-neg-file-eval is given and --keep-ann is False, only evaluate GO terms that have at least 2% of annotations. Useful to speed-up processing for term-based algorithms")
    parser.add_option_group(group)

    run_algs.add_alg_opts(parser)
    run_algs.add_string_opts(parser)

    # additional parameters
    group = OptionGroup(parser, 'Additional options')
    group.add_option('-T', '--taxon', type='string', action='append',
                      help="Specify the species taxonomy ID to use to be left out. Multiple may be specified. Otherwise, all species will be used")
    group.add_option('', '--postfix', type='string', default="",
                     help="Postfix to append to output file. Useful if running multiple in parallel. TODO figure out how to automatically combine the multiple output files")  # for file to write. Otherwise the standard outputs/<version>/all/<algorithm>/<exp_name>/ will be used. 
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
        opts.version = ["2018_06-seq-sim-1e-25"]

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
    if 'aptrank' in opts.algorithm:
        from aptrank.aptrank import AptRank

    if opts.keep_ann and opts.pos_neg_file_eval is None:
        print("ERROR: Must specify a --pos-neg-file-eval to use the --keep-ann option")
        sys.exit(1)

    opts.k = opts.k if opts.k is not None else [200]
    opts.t = opts.t if opts.t is not None else [2]
    opts.s = opts.s if opts.s is not None else [200]
    # default for deltaUBLB is None
    opts.eps = opts.eps if opts.eps is not None else [0.0001]
    opts.epsUB = opts.epsUB if opts.epsUB is not None else [0]

    run_algs.validate_string_opts(opts)

    return opts


def main(version, exp_name, W, prots, ann_matrix, goids,
         algorithms, opts, taxons=None, eval_ann_matrix=None):
    """
    *W*: If using STRING networks, then W should be a tuple containing
        a list of all STRING sparse matrices to use when weighting (SWSN)
        a list of the names of the STRING networks
    *eval_ann_matrix*: matrix with the same data as *ann_matrix*, except 
        it will only be used for evaluation
    """
    # option to use negative examples when evaluating predictions
    #use_negatives_for_eval = True 
    # start with less GO terms
    #exp_name = "eval-species-%s-%d-%d-25iters" % (evidence_codes, cut1, cut2)
    if opts.non_pos_as_neg_eval is False:
        exp_name += "-use-neg" 
    else:
        exp_name += "-non-pos-neg" 
    if opts.oracle:
        exp_name += "-oracle" 
    if opts.keep_ann:
        exp_name += "-keep-ann" 

    # TODO make this more streamlined
    aptrank_data = None 
    if 'birgrank' in algorithms or 'aptrank' in algorithms:
        dag_matrix, pos_matrix, dag_goids = run_birgrank.setup_h_ann_matrices(
                prots, opts.obo_file, opts.pos_neg_file, goterms=goids)
        aptrank_data = (dag_matrix, pos_matrix, dag_goids)

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
    global species_to_uniprot_idx
    species_to_uniprot_idx = defaultdict(set)
    for p in uniprot_to_species:
        species_to_uniprot_idx[uniprot_to_species[p]].add(node2idx.get(p))
    for s in species_to_uniprot_idx:
        species_to_uniprot_idx[s].discard(None) 

    print("Training/testing with %d species, %d goterms" % (len(taxons), len(goids)))
    # for each alg, taxon and go term pair, see which already exist and skip them
    alg_taxon_terms_to_skip = {alg: {} for alg in algorithms}

    for alg in algorithms:
        #out_file = "%s/all/%s/%s/loso-%s%s%sl%d-a%s-eps%s-maxi%d%s%s%s%s.txt" % (  # for now, look at the main file to see which are missing
        out_file = "%s/all/%s/%s/loso-%s%s%sl%d-a%s-eps%s-maxi%d%s%s%s.txt" % (
            RESULTSPREFIX, alg, exp_name,
            'unw-' if opts.unweighted else '', 
            'goterm-weight-' if opts.weight_per_goterm else '',
            'swsn-' if opts.weight_swsn else '',
            0 if opts.sinksourceplus_lambda is None else opts.sinksourceplus_lambda,
            str(opts.alpha).replace('.', '_'),
            str(opts.eps).replace('.', '_'),
            opts.max_iters,
            '-t%s-m%s-l%s' % (
                str(opts.theta).replace('.','_'), str(opts.mu).replace('.','_'),
                str(opts.br_lambda).replace('.','_')) if alg == 'birgrank' else '',
            '-l%s-k%s-s%s-t%s-%s' % (
                str(opts.br_lambda).replace('.','_'), str(opts.apt_k).replace('.', '_'),
                str(opts.apt_s).replace('.', '_'), str(opts.apt_t).replace('.', '_'),
                opts.diff_type) if alg == 'aptrank' else '',
            '-tol%s' % (str(opts.tol).replace('.', '_')) if alg == 'genemania' else '',
            #opts.postfix,
        )

        if os.path.isfile(out_file) and opts.forcealg:
            if len(taxons) > 1: 
                print("Removing %s as results will be appended to it for each taxon" % (out_file))
                os.remove(out_file)
        elif os.path.isfile(out_file) and opts.only_pred is False:
            if len(taxons) > 1: 
                print("WARNING: %s results file already exists. Appending to it" % (out_file))
                #print("%s results file already exists. Use --forcealg to overwrite it.\nQuitting" % (out_file))
                #sys.exit()
                # check which results already exist and append to the rest
                print("Reading results from %s " % (out_file))
                taxon_terms_completed = utils.readColumns(out_file, 1, 2)
                alg_taxon_terms_to_skip[alg] = {taxon: set() for taxon, term in taxon_terms_completed}
                for taxon, term in taxon_terms_completed:
                    alg_taxon_terms_to_skip[alg][taxon].add(term)
                print("\t%d taxon - term pairs already finished" % (len(taxon_terms_completed)))
        elif opts.only_pred is False:
            print("Writing results to %s" % (out_file))

    # convert the annotation matrix to a lil matrix
    # Update: This uses much more RAM
    #ann_matrix = ann_matrix.tolil()
    #if eval_ann_matrix is not None:
    #    eval_ann_matrix = eval_ann_matrix.tolil()
    goid_scores = None
    if opts.keep_ann:
        print("\nRunning algs with all annotations in the --pos-neg-file and evaluating based on the annotations in the --pos-neg-file-eval.")
        print("\tindividual species will be evaluated after")
        train_ann_mat, test_ann_mat, sp_goterms = leave_out_taxon(
            None, ann_matrix, goids, prots, species_to_uniprot_idx,
            eval_ann_matrix=eval_ann_matrix, keep_ann=opts.keep_ann,
            non_pos_as_neg_eval=opts.non_pos_as_neg_eval, 
            oracle=opts.oracle)
        print("\t%d/%d goterms have at least 1 train_pos and train_neg" % (len(sp_goterms), len(goids)))
        goid_scores, params_results = run_and_eval_algs(version, exp_name, W, prots, goids,
                          algorithms, opts, train_ann_mat, test_ann_mat,
                          taxons=taxons, aptrank_data=aptrank_data, taxon='all')

    else:
        params_results = defaultdict(int)
        # split the sets of positives/negatives into 19 sets with one species left out of each
        for s in tqdm(sorted(taxons)):
            tqdm.write("\n" + "-"*30)
            #print("\n" + "-"*30)
            #print("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
            tqdm.write("Taxon: %s - %s" % (
                s, selected_species[s]))

            train_ann_mat, test_ann_mat, sp_goterms = leave_out_taxon(
                s, ann_matrix, goids, prots, species_to_uniprot_idx,
                eval_ann_matrix=eval_ann_matrix, keep_ann=opts.keep_ann,
                non_pos_as_neg_eval=opts.non_pos_as_neg_eval, 
                eval_goterms_with_left_out_only=opts.eval_goterms_with_left_out_only,
                oracle=opts.oracle, cutoff=opts.num_test_cutoff,
                terms_to_skip=alg_taxon_terms_to_skip[alg][s] if s in alg_taxon_terms_to_skip[alg] else None)

            tqdm.write("\t%d/%d goterms with >= %d annotations" % (len(sp_goterms), len(goids), opts.num_test_cutoff))
            if len(sp_goterms) == 0:
                print("\tskipping")
                continue

            goid_scores, curr_params_results = run_and_eval_algs(
                version, exp_name, W, prots, goids,
                algorithms, opts, train_ann_mat, test_ann_mat,
                taxons=taxons, aptrank_data=aptrank_data, taxon=s)
            # don't store the scores for all taxons. Takes too much space
            for key in curr_params_results:
                params_results[key] += curr_params_results[key]

    print("Final running times:")
    print(", ".join(["%s: %0.4f" % (key, val) for key, val in sorted(params_results.items())]))
    # return the goid_scores even though they will only be for a single taxon
    # because other scripts may use 
    return goid_scores, params_results


def run_and_eval_algs(
        version, exp_name, W, prots, goids,
        algorithms, opts, train_ann_mat, test_ann_mat,
        taxons=None, aptrank_data=None, taxon="all"):
    params_results = defaultdict(int)

    if 'birgrank' in algorithms or 'aptrank' in algorithms:
        dag_matrix, pos_matrix, dag_goids = aptrank_data 
        # TODO the matrix for BirgRank and the train_matrix do not have the same goids. 
        # I need to build a pos_mat with the train_matrix annotations
        #train_pos_mat = sparse.lil_matrix(pos_matrix.shape)
        # TODO test this
        train_pos_mat = sparse.csr_matrix(pos_matrix.shape)
        dag_goids2idx = {g: i for i, g in enumerate(dag_goids)}
        for i in range(len(goids)):
            dag_goid_idx = dag_goids2idx[goids[i]]
            train_pos_mat[dag_goid_idx] = train_ann_mat[i]
        # now set the negatives to 0 as birgrank doesn't use negatives
        train_pos_mat[train_pos_mat < 0] = 0
        # not needed for lil matrix
        #train_pos_mat.eliminate_zeros() 
        # Also, speed-up birgrank by only getting the scores for the nodes that will be evaluated (a positive or negative for at least 1 GO term)
        test_nodes = set()
        for i in range(len(goids)):
            test_nodes.update(set(list(test_ann_mat[i].nonzero()[1])))

    #if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
    if opts.weight_swsn:
        # use the simultaneous weighting method to weight the networks
        sparse_networks, net_names = W
        # commented out the out_file as there is no need to keep storing them
        #out_file = "inputs/%s/%s/%d-nets-combined-SWSN-leave-out-%s.npz" % (
        #    version, exp_name, len(sparse_networks), s)
        #if os.path.isfile(out_file):
        #    print("Loading SWSN weighted network from %s" % (out_file))
        #    W = sparse.load_npz(out_file)
        #else:
        W, swsn_time = setup.weight_SWSN(train_ann_mat, sparse_networks,
                net_names=net_names, nodes=prots)
        print("\ttotal time to weight networks swsn: %0.4f sec" % (swsn_time))
        params_results['swsn_time'] += swsn_time
                #net_names=net_names, out_file=out_file, nodes=prots)

    if opts.keep_ann:
        print("Keeping all annotations when making predictions")
    if opts.only_pred:
        print("Making predictions only and writing %d to a file" % (opts.num_pred_to_write))
        test_ann_mat = None
    elif opts.non_pos_as_neg_eval is True: 
        print("Evaluating using all non-ground-truth positives for the taxon as false positives")
    else:
        print("Evaluating using only the ground-truth negatives predicted as positives as false positives")

    for alg in algorithms:
        rank_pos_neg = test_ann_mat if opts.rank_pos_neg is True else None
        # for now, use most of the defaults
        # change to running sinksource with 25 iterations
        # leave alpha at default 0.8
        alg_runner = run_algs.Alg_Runner(
            opts.eng, version, exp_name, W, prots, train_ann_mat, goids,
            algorithms=[alg], weight_swsn=False, weight_per_goterm=opts.weight_per_goterm, unweighted=opts.unweighted,
            ss_lambda=opts.sinksourceplus_lambda, taxon=taxon,
            eps=opts.eps, epsUB_list=opts.epsUB, max_iters=opts.max_iters,
            rank_topk=opts.rank_topk, rank_all=opts.rank_all, rank_pos_neg=rank_pos_neg, compare_ranks=opts.compare_ranks,  # sinksource-squeeze parameters
            tol=opts.tol,  # genemania parameters
            k_list=opts.k, t_list=opts.t, s_list=opts.s,  # sinksource-ripple parameters
            alpha=opts.alpha, theta=opts.theta, mu=opts.mu, br_lambda=opts.br_lambda,  # birgrank parameters
            k=opts.apt_k, s=opts.apt_s, t=opts.apt_t, diff_type=opts.diff_type,  # aptrank parameters
            num_pred_to_write=opts.num_pred_to_write, verbose=opts.verbose, 
            forcealg=opts.forcealg, progress_bar=False)
        if alg == 'birgrank' or alg == 'aptrank':
            # the W matrix is already normalized, so I can run birgrank/aptrank from here
            goid_scores, curr_params_results = alg_runner.run_aptrank_with_params(
                train_pos_mat, dag_matrix, alg=alg, nodes=test_nodes) 
            curr_goids = list(dag_goids)
        else:
            curr_goids = list(goids)  # .copy() doesn't work for some reason
            goid_scores, curr_params_results = alg_runner.main()
        #print(curr_params_results)
        for key in curr_params_results:
            params_results[key] += curr_params_results[key]

        # skip evaluation if --only-pred is specified
        if opts.only_pred:
            continue
        # now evaluate 
        out_dir = "outputs/%s/all/%s/%s/" % (version, alg, exp_name)

        # this will write an file containing the fmax for each goterm 
        # with the taxon name in the name of the file
        write_prec_rec = False 
        if len(taxons) == 1:
            print("Also writing prec/rec stats")
            write_prec_rec = True 
            out_dir += "goids/"

        utils.checkDir(out_dir)
        out_pref = "%s/%sloso-%s%s%sl%d-a%s-eps%s-maxi%d" % (
            out_dir, "all-sp-" if taxon=='all' else '',
            'unw-' if opts.unweighted else '',
            'goterm-weight-' if opts.weight_per_goterm else '',
            'swsn-' if opts.weight_swsn else '',
            0 if alg_runner.ss_lambda is None else int(alg_runner.ss_lambda),
            str(opts.alpha).replace('.','_'), str(opts.eps).replace('.','_'),
            opts.max_iters)
        if alg == 'birgrank':
            out_pref += '-t%s-m%s-l%s' % (
                str(opts.theta).replace('.','_'), str(opts.mu).replace('.','_'), str(opts.br_lambda).replace('.','_'))
        elif alg == 'aptrank':
            k, s, t, diff_type = opts.apt_k, opts.apt_s, opts.apt_t, opts.diff_type
            out_pref += 'l%s-k%s-s%s-t%s-%s' % (
                str(opts.br_lambda).replace('.', '_'), str(k).replace('.', '_'), str(s).replace('.', '_'),
                str(t).replace('.', '_'), diff_type)
        elif alg == 'genemania':
            out_pref += '-tol%s' % (str(opts.tol).replace('.', '_') if alg == 'genemania' else '')

        out_file = "%s%s.txt" % (out_pref, opts.postfix)
        alg_runner.evaluate_ground_truth(
            goid_scores, curr_goids, test_ann_mat, out_file,
            #non_pos_as_neg_eval=opts.non_pos_as_neg_eval,
            taxon=taxon, write_prec_rec=write_prec_rec, 
            #append=True if taxon != 'all' else False)  # update: to parallelize keep_ann, always append
            append=True)

        #if opts.keep_ann and per_taxon is True:
        if opts.keep_ann:
            out_pref += "-per-taxon"
            out_file = "%s%s.txt" % (out_pref, opts.postfix)
            # now get the fmax per taxon
            print("Getting the per-taxon fmax results")
            all_train_ann_mat, all_test_ann_mat = train_ann_mat, test_ann_mat
            for s in taxons:
                train_ann_mat, test_ann_mat, sp_goterms = leave_out_taxon(
                    s, all_train_ann_mat, goids, prots, species_to_uniprot_idx,
                    eval_ann_matrix=all_test_ann_mat, keep_ann=opts.keep_ann,
                    non_pos_as_neg_eval=opts.non_pos_as_neg_eval, 
                    oracle=opts.oracle, cutoff=opts.num_test_cutoff)

                tqdm.write("\t%d/%d goterms with >= %d annotations" % (len(sp_goterms), len(goids), opts.num_test_cutoff))
                if len(sp_goterms) == 0:
                    print("\tskipping")
                    continue
                alg_runner.evaluate_ground_truth(
                    goid_scores, curr_goids, test_ann_mat, out_file,
                    #non_pos_as_neg_eval=opts.non_pos_as_neg_eval,
                    taxon=s, write_prec_rec=write_prec_rec, 
                    append=True)

    return goid_scores, params_results


def leave_out_taxon(s, ann_matrix, goids, prots, species_to_uniprot_idx,
                    eval_ann_matrix=None, keep_ann=False, 
                    non_pos_as_neg_eval=False, eval_goterms_with_left_out_only=False,
                    oracle=False, cutoff=1, terms_to_skip=None):
    """
    Training positives are removed from testing positives, and train pos and neg are removed from test neg
        I don't remove training negatives from testing positives, because not all algorithms use negatives
    *s*: species to be left out. If s is None, then no species will be left out, and keep_ann must be True.
    *ann_matrix*: should be a lil matrix
    *eval_ann_matrix*: should be a lil matrix
    *eval_goterms_with_left_out_only*: if eval_ann_matrix is given and keep_ann is False, 
        only evaluate GO terms that have at least 2% of annotations. 
        Useful to speed-up processing for term-based algorithms
    *oracle*: remove train negatives that are actually test positives
    *cutoff*: minimum number of annotations for each go term in the left-out species 
    *alg_taxon_terms_to_skip*: set of terms to be skipped 
    """
    # leave this taxon out by removing its annotations
    # rather than a dictionary, build a matrix
    train_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
    test_ann_mat = sparse.lil_matrix(ann_matrix.shape, dtype=np.float)
    sp_goterms = []
    skipped_eval_no_left_out_ann = 0
    #pdb.set_trace()
    for i in range(len(goids)):
        if terms_to_skip is not None and goids[i] in terms_to_skip:
            continue
        pos, neg = alg_utils.get_goid_pos_neg(ann_matrix, i)
        ann_pos = set(list(pos))
        ann_neg = set(list(neg))
        # first setup the training annotations (those used as positives/negatives for the algorithm)
        if keep_ann:
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
            # if this species has little-to-no annotations that are being left-out, then we can skip it
            if not keep_ann and eval_goterms_with_left_out_only:
                # make sure we don't divide by zero
                if len(train_pos) < cutoff:
                    continue
                # If the percentage of left-out ann is less than 1%, then skip it
                if (len(ann_pos) - len(train_pos)) / float(len(train_pos)) < .02:
                    skipped_eval_no_left_out_ann += 1 
                    continue
                #else:
                #    print(len(ann_pos), len(train_pos), (len(ann_pos) - len(train_pos)) / float(len(train_pos)))
        if s is None:
            test_pos = eval_pos
            test_neg = eval_neg
            if non_pos_as_neg_eval:
                # everything minus the positives
                test_neg = set(prots) - test_pos
        else:
            # TODO I should limit these to the proteins in the network
            test_pos = eval_pos & species_to_uniprot_idx[s]
            # UPDATE 2018-06-27: Only evaluate the species prots as negatives, not all prots
            if non_pos_as_neg_eval:
                test_neg = species_to_uniprot_idx[s] - eval_pos
                test_neg.discard(None)
                #test_neg = np.asarray(sorted(test_neg)).astype(int)
            else:
                test_neg = eval_neg & species_to_uniprot_idx[s]
        # UPDATE 2018-06-30: Remove test positives/negatives that are part of the training positives/negatives
        # don't remove test positives if its a training negative because not all algorithms use negatives
        test_pos -= train_pos 
        if oracle:
            train_neg -= test_pos
        test_neg -= train_pos | train_neg 
        # UPDATE 2018-10: Add a cutoff on both the # of training positive and # of test pos
        if len(train_pos) < cutoff or len(test_pos) < cutoff or \
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

    if eval_ann_matrix is not None and not keep_ann and eval_goterms_with_left_out_only:
        print("\t%d goterms skipped_eval_no_left_out_ann (< 0.02 train ann in the left-out species)" % (skipped_eval_no_left_out_ann))

    return train_ann_mat.tocsr(), test_ann_mat.tocsr(), sp_goterms


def run(args):
    opts = parse_args(args)

    opts.eng = run_algs.get_matlab_engine()
    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 

    #if 'birgrank' in opts.algorithm and len(opts.pos_neg_file) > 1:
    #    print("Birgrank not yet implemented with multiple hierarchies. Use only bp or mf")
    #    sys.exit()

    #goid_pos, goid_neg = alg_utils.parse_pos_neg_files(opts.pos_neg_file, goterms=goterms) 
    # load the network matrix and protein IDs
    #net_file = opts.net_file
    #if net_file is None:
    INPUTSPREFIX, _, net_file, selected_strains = f_settings.set_version(opts.version) 
    # TODO this should be better organized so that any STRING networks
    # can be used
    #if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
    if opts.weight_swsn or opts.weight_per_goterm:
        out_pref_net = "%s/sparse-nets/" % (INPUTSPREFIX)
        utils.checkDir(out_pref_net)
        # build the file containing the sparse networks
        sparse_networks, network_names, prots = setup.create_sparse_net_file(
            opts.version, out_pref_net, selected_strains=selected_strains,
            string_nets=opts.string_networks, string_cutoff=f_settings.VERSION_STRING_CUTOFF[opts.version],
            forcenet=False)
        # TODO organize this better
        W = (sparse_networks, network_names)
    else:
        W, prots = alg_utils.setup_sparse_network(net_file)

    # now build the annotation matrix
    ann_matrix, goids = setup.setup_sparse_annotations(opts.pos_neg_file, goterms, prots)
    #print(goids)
    test_ann_matrix = None 
    if opts.pos_neg_file_eval is not None:
        test_ann_matrix, test_goids = setup.setup_sparse_annotations(opts.pos_neg_file_eval, set(goids), prots)
        # Make sure the goids are the same before using it
        num_ann_with_test = len(set(goids) & set(test_goids))
        if num_ann_with_test != 0:
            #print(test_goids)
            print("WARNING: %d / %d goids in the annotation matrix are in the evaluation matrix" % (num_ann_with_test, len(goids)))
            ann_goids_without_eval = set(goids) - set(test_goids)

            print("\tremoving %d GO terms from the annotation matrix" % (len(ann_goids_without_eval)))
            goids_to_remove = []
            new_goids = []
            for i, g in enumerate(goids):
                if g in ann_goids_without_eval:
                    goids_to_remove.append(i)
                    continue
                new_goids.append(g)
            ann_matrix = alg_utils.delete_rows_csr(ann_matrix.tocsr(), goids_to_remove)
            goids = new_goids

        print("\tmatching the evaluation matrix to the annotation matrix")
        test_goids2idx = {g: i for i, g in enumerate(test_goids)}
        # slim down the evaluation matrix to have the same GO terms as the ann matrix
        #matching_test_matrix = sparse.csr_matrix(ann_matrix.shape)
        matching_test_matrix = sparse.lil_matrix(ann_matrix.shape)
        for i, g in enumerate(tqdm(goids)):
            matching_test_matrix[i] = test_ann_matrix[test_goids2idx[g]]
        test_ann_matrix = matching_test_matrix.tocsr()
        del matching_test_matrix
    return opts, W, prots, ann_matrix, goids, test_ann_matrix


if __name__ == "__main__":
    opts, W, prots, ann_matrix, goids, test_ann_matrix = run(sys.argv)
    # TODO alpha doesn't change for all algorithms, just SS
    main(opts.version, opts.exp_name, W, prots, ann_matrix, goids,
         opts.algorithm, opts, taxons=opts.taxon,
         eval_ann_matrix=test_ann_matrix)

