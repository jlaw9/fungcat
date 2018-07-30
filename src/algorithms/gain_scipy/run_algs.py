
# Quick script to run/test the algorithms
#print("Importing libraries")

from optparse import OptionParser,OptionGroup
from collections import defaultdict
import os
import sys
from tqdm import tqdm
import itertools
sys.path.append("src")
import utils.file_utils as utils
# add the folder above and two folders above
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fungcat_settings as f_settings
import setup_sparse_networks as setup
import alg_utils
import sinksource
import genemania
#import sinksource_topk_ub
import sinksource_ripple
import sinksource_squeeze
from aptrank.birgrank import birgRank
import aptrank.run_birgrank as run_birgrank
#import pandas as pd
import networkx as nx
from scipy import sparse
import numpy as np
from sklearn import metrics
from sklearn.model_selection import KFold
#import gc


ALGORITHMS = [
    "sinksourceplus-ripple",  # This uses the same UB as Ripple
    "sinksource-ripple",  # This uses the same UB as Ripple, but with negatives
#    "sinksource-topk-ub",  # This computes UB and LB using iteration
#    "sinksourceplus-topk-ub",  # This uses the same UB as sinksource-ub-topk
    "sinksourceplus-squeeze",
    "sinksource-squeeze",
    "sinksourceplus",  # same as sinksource-ovn
    "sinksource",  # same as sinksource-ova, but with the option of using lambda and/or alpha
    "localplus",  # same as local-ovn, but with the option of using lambda and/or alpha
    "local",  # same as local-ova, but with the option of using lambda and/or alpha
    "birgrank",
    "genemania",
    ]

# These algorithms don't use any negative examples
POSITIVE_ALGS = [
    "sinksourceplus-ripple",  
    "sinksourceplus-squeeze",
    "sinksourceplus",  
    "localplus",  
    ]

class Alg_Runner:
    """
    Base class for running algorithms
    """
    def __init__(
            self, version, exp_name,
            W, prots, ann_matrix, goids,
            algorithms=["sinksource"], weight_swsn=False,
            unweighted=False, ss_lambda=None, 
            alpha=0.8, eps=0.0001, max_iters=1000,
            k_list=[100], t_list=[2], s_list=[50],  
            deltaUBLB_list=[None], epsUB_list=[0], 
            rank_topk=False, rank_all=False, rank_pos=False, compare_ranks=False,
            taxon=None, num_pred_to_write=100, 
            aptrank_data=None, theta=.5, mu=.5,
            only_cv=False, cross_validation_folds=None, 
            forcenet=False, forcealg=False, verbose=False, progress_bar=True):
        """
        *eps*: Convergence cutoff for sinksource and sinksourceplus
        *aptrank_data*: tuple containing the dag_matrix, annotation matrix, and the goids of the hierarchy.
        """
        self.version = version
        self.exp_name = exp_name
        #self.W = W
        self.prots = prots
        self.ann_matrix = ann_matrix.tolil()
        self.goids = goids
        if aptrank_data is not None:
            self.dag_matrix, self.pos_matrix, self.dag_goids = aptrank_data
        self.theta, self.mu = theta, mu
        self.algorithms = algorithms
        self.weight_swsn = weight_swsn
        self.unweighted = unweighted
        self.ss_lambda = ss_lambda
        # parameters for algorithms
        self.alpha, self.eps, self.max_iters = alpha, eps, max_iters
        # ripple/squeeze parameters
        self.k_list = k_list
        self.t_list = t_list
        self.s_list = s_list
        self.deltaUBLB_list = deltaUBLB_list
        self.epsUB_list = epsUB_list
        self.rank_topk = rank_topk
        self.rank_all = rank_all
        self.rank_pos = rank_pos
        self.compare_ranks = compare_ranks
        # evaluation parameters
        self.taxon = taxon
        self.only_cv = only_cv
        self.cross_validation_folds = cross_validation_folds
        self.num_pred_to_write = num_pred_to_write
        self.forcenet = forcenet
        self.forcealg = forcealg
        self.verbose = verbose
        # option to display progress bar while iterating over GO terms
        self.progress_bar = progress_bar

        if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[self.version] and not self.unweighted:
            print("\tkeeping individual STRING networks to be weighted later")
            self.sparse_networks, self.net_names = W
            print("\tnormalizing the networks")
            self.normalized_nets = []
            for net in self.sparse_networks:
                self.normalized_nets.append(setup._net_normalize(net))
        else:
            if self.unweighted is True:
                print("\tsetting all edge weights to 1 (unweighted)")
                #and re-normalizing by dividing each edge by the node's degree")
                W = (W > 0).astype(int) 
            self.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=self.ss_lambda)

            if 'genemania' in self.algorithms:
                print("\nCreating Laplacian matrix for GeneMANIA")
                self.L = genemania.setup_laplacian(W)

        # used to map from node/prot to the index and vice versa
        #self.idx2node = {i: n for i, n in enumerate(prots)}
        self.node2idx = {n: i for i, n in enumerate(prots)}
        # used to map from index to goid and vice versa
        #self.idx2goid = {i: g for i, g in enumerate(goid)}
        self.goid2idx = {g: i for i, g in enumerate(goids)}

    def main(self):
        INPUTSPREFIX, RESULTSPREFIX, _, selected_strains = f_settings.set_version(self.version)
        self.INPUTSPREFIX = INPUTSPREFIX
        self.RESULTSPREFIX = RESULTSPREFIX
        # start with a single species
        #self.taxon = "208964"
        #if self.taxon is not None:
        #    network_file = f_settings.STRING_TAXON_UNIPROT % (self.taxon, self.taxon, f_settings.STRING_CUTOFF)
        #    gain_file = f_settings.FUN_FILE % (self.taxon, self.taxon)
        # organize the results by algorithm, then algorithm parameters, then GO term
        #params_results = {}
        for alg in self.algorithms:
            self.out_dir = "%s/all/%s/%s" % (self.RESULTSPREFIX, alg, self.exp_name)
            utils.checkDir(self.out_dir)

            # TODO add an option to run birgrank in prediction mode
            if alg == 'birgrank' and self.only_cv is False:
                continue

            out_pref = "%s/pred-%s%s" % (
                self.out_dir, 'unw-' if self.unweighted else '',
                'l'+str(self.ss_lambda).replace('.', '_')+'-' if self.ss_lambda is not None else '')
            # "prediction mode" is run if the user doesn't specify --only-cv
            if self.only_cv is False:
                if self.num_pred_to_write == 0:
                    out_pref = None

                # now run the algorithm with all combinations of parameters
                # TODO the goid scores are for a single set of parameters. 
                # TODO the goid scores are only stored if the ground_truth is given
                # as storing all goid_scores uses too much RAM
                goid_scores, curr_params_results = self.run_alg_with_params(
                        alg, out_pref=out_pref)

                #params_results.update(curr_params_results)

            # "cross validation mode" is run if the user specifies --cross-validation-folds X
            if self.cross_validation_folds is not None:
                if self.weight_swsn is True or alg == "birgrank":
                    self.run_cv_all_goterms(alg, folds=self.cross_validation_folds)
                else:
                    out_pref = "%s/cv-%dfolds-%sl%d-" % (self.out_dir, self.cross_validation_folds,
                        'unw-' if self.unweighted else '', 0 if self.ss_lambda is None else int(self.ss_lambda))
                    cv_params_results = self.run_goterm_cv(alg,
                            folds=self.cross_validation_folds, out_pref=out_pref)
                    #params_results.update(cv_params_results)

            #if self.verbose:
            #    print(params_results_to_table(params_results))

        #version_params_results[version] = params_results

        # I implemented this for Ripple and Squeeze
        # TODO move this somewhere else
#        #for version in self.versions:
#        table = params_results_to_table(params_results)
#        if self.verbose:
#            print(self.version)
#            print(table)
#        out_dir = "outputs/viz/params-results"
#        utils.checkDir(out_dir)
#        out_file = "%s/%s-%s-params-results.tsv" % (out_dir, self.exp_name, self.version)
#        if self.forcealg is True or not os.path.isfile(out_file):
#            print("Writing params-results to: %s" % (out_file))
#            with open(out_file, 'w') as out:
#                out.write(table)
#        else:
#            print("%s already exists. Use --forcealg to overwrite it" % (out_file))

        if self.only_cv is False:
            return goid_scores
        else:
            return

    def run_cv_all_goterms(self, alg, folds=5):
        """
        Split the positives and negatives into folds across all GO terms
        and then run the algorithms on those folds.
        For now, only setup to run BirgRank. TODO allow running each algorithm on the same split of data. 
        """
        # because each fold contains a different set of positives, and combined they contain all positives,
        # store all of the prediction scores from each fold in a matrix
        combined_fold_scores = sparse.lil_matrix(self.ann_matrix.shape, dtype=np.float)
        orig_ann_matrix = self.ann_matrix
        out_pref = "%s/cv-%dfolds-all-goids-%sl%d-" % (self.out_dir, self.cross_validation_folds,
                'unw-' if self.unweighted else '', 0 if self.ss_lambda is None else int(self.ss_lambda))
        if alg == 'birgrank':
            alpha, theta, mu = self.alpha, self.theta, self.mu
            out_pref += 'a%s-t%s-m%s' % (
                str(alpha).replace('.','_'), str(theta).replace('.','_'),
                str(mu).replace('.','_'))
            dag_goids2idx = {g: i for i, g in enumerate(self.dag_goids)}
        utils.checkDir(os.path.dirname(out_pref))

        ann_matrix_folds = self.split_cv_all_goterms(folds=folds)
        for curr_fold, (train_ann_mat, test_ann_mat) in enumerate(ann_matrix_folds):
            print("*  "*20)
            print("Fold %d" % (curr_fold+1))

            # weight the network if needed
            if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[self.version] and not self.unweighted:
                W = setup.weight_SWSN(train_ann_mat, self.sparse_networks,
                        net_names=self.net_names, nodes=self.prots)
                self.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=self.ss_lambda)
                if alg == 'genemania':
                    self.L = genemania.setup_laplacian(W)

            # run each algorithm on these folds
            if alg == 'birgrank':
                # TODO the matrix for BirgRank and the train_matrix do not have the same goids. 
                # I need to build a pos_mat with the train_matrix annotations
                train_pos_mat = sparse.lil_matrix(self.pos_matrix.shape)
                # the birgrank/aptrank annotation (positive) matrix has a different # GO terms,
                # but the same node columns. So I just need to align the rows
                # Also, speed-up birgrank by only getting the scores for the nodes that are a positive or negative for at least 1 GO term
                test_nodes = set()
                for i in range(len(self.goids)):
                    dag_goid_idx = dag_goids2idx[self.goids[i]]
                    train_pos_mat[dag_goid_idx] = train_ann_mat[i]
                    test_nodes.update(set(list(train_ann_mat[i].nonzero()[1])))
                # now set the negatives to 0 as birgrank doesn't use negatives
                train_pos_mat[train_pos_mat < 0] = 0

                # the W matrix is already normalized, so I can run
                # birgrank/aptrank from here
                Xh = self.run_aptrank_with_params(
                    train_pos_mat, self.dag_matrix, alg=alg, nodes=test_nodes) 
                goid_scores = sparse.lil_matrix(self.ann_matrix.shape, dtype=np.float)
                # limit the scores to only the GOIDs for which we have annotations
                for i in range(len(self.goids)):
                    goid_scores[i] = Xh[dag_goids2idx[self.goids[i]]]
            else:
                # change the annotation matrix to the current fold
                self.ann_matrix = train_ann_mat
                goid_scores, _ = self.run_alg_with_params(alg)

            # store only the scores of the test (left out) positives and negatives
            for i in range(len(self.goids)):
                test_pos, test_neg = alg_utils.get_goid_pos_neg(test_ann_mat, i)
                curr_goid_scores = goid_scores[i].toarray().flatten()
                curr_comb_scores = combined_fold_scores[i].toarray().flatten()
                curr_comb_scores[test_pos] = curr_goid_scores[test_pos]
                curr_comb_scores[test_neg] = curr_goid_scores[test_neg]
                combined_fold_scores[i] = curr_comb_scores 

        #curr_goids = self.dag_goids if alg == 'birgrank' else self.goids
        # now evaluate the results and write to a file
        self.evaluate_ground_truth(
            combined_fold_scores, self.goids, orig_ann_matrix, out_pref,
            #non_pos_as_neg_eval=opts.non_pos_as_neg_eval,
            alg=alg)

    def split_cv_all_goterms(self, folds=5):
        """
        Split the positives and negatives into folds across all GO terms
        Required for running CV with birgrank/aptrank
        *returns*: a list of tuples containing the (train pos, train neg, test pos, test neg)
        """
        print("Splitting all annotations into %d folds by splitting each GO terms annotations into folds, and then combining them" % (folds))
        # TODO there must be a better way to do this than getting the folds in each go term separately
        # list of tuples containing the (train pos, train neg, test pos, test neg) 
        ann_matrix_folds = []
        for i in range(folds):
            train_ann_mat = sparse.lil_matrix(self.ann_matrix.shape, dtype=np.float)
            test_ann_mat = sparse.lil_matrix(self.ann_matrix.shape, dtype=np.float)
            ann_matrix_folds.append((train_ann_mat, test_ann_mat))

        for i in tqdm(range(self.ann_matrix.shape[0]), total=self.ann_matrix.shape[0], disable=not self.progress_bar):
            goid = self.goids[i]
            positives, negatives = alg_utils.get_goid_pos_neg(self.ann_matrix, i)
            #print("%d positives, %d negatives for goterm %s" % (len(positives), len(negatives), goid))
            kf = KFold(n_splits=folds, shuffle=True)
            kf_neg = KFold(n_splits=folds, shuffle=True)
            kf.get_n_splits(positives)
            kf_neg.get_n_splits(negatives)
            fold = 0

            for (pos_train_idx, pos_test_idx), (neg_train_idx, neg_test_idx) in zip(kf.split(positives), kf_neg.split(negatives)):
                train_pos, test_pos= positives[pos_train_idx], positives[pos_test_idx]
                train_neg, test_neg = negatives[neg_train_idx], negatives[neg_test_idx]
                train_ann_mat, test_ann_mat = ann_matrix_folds[fold]
                fold += 1
                # build an array of the scores and set it in the goid sparse matrix of scores
                for pos, neg, mat in [(train_pos, train_neg, train_ann_mat), (test_pos, test_neg, test_ann_mat)]:
                    pos_neg_arr = np.zeros(len(self.prots))
                    pos_neg_arr[list(pos)] = 1
                    pos_neg_arr[list(neg)] = -1
                    mat[i] = pos_neg_arr

        return ann_matrix_folds

    def run_aptrank_with_params(self, pos_mat, hierarchy_mat, alg='birgrank', nodes=None, out_pref=None):
        """ Run a protein-based algorithm that uses the hierarchy
        *hierarchy_matrix*: matrix of hierarchy relationships
        *pos_matrix*: matrix of goterm - protein annotations (1 or 0).
            Normally these would be the leaf annotations, not propagated
        *nodes*: set of nodes for which to run RWR and get GO term scores

        *returns*: a matrix of prediction scores with the same
            dimensions as pos_mat
        """

        assert (pos_mat.shape[0] == hierarchy_mat.shape[0]), \
            "Error: annotation and hierarchy matrices " + \
            "do not have the same shape: %d, %d" % (
                pos_mat.shape[0], hierarchy_mat.shape[0])

        # remove the negatives as aptrank doesn't use them
        #pos_mat = self.ann_matrix.copy()
        #pos_mat[pos_mat < 0] = 0
        if alg == 'birgrank':
            Xh = birgRank(self.P, pos_mat.transpose(), hierarchy_mat,
                        alpha=self.alpha, theta=self.theta, mu=self.mu, 
                        eps=self.eps, max_iters=self.max_iters,
                        nodes=nodes, verbose=self.verbose)
            Xh = Xh.T
        else:
            print("alg %s not yet implemented." % (alg))
            return 
        # now write the scores to a file
        if out_pref is not None:
            out_file = "%sa%s-t%s-m%s.txt" % (
                out_pref, str(self.alpha).replace('.', '_'),
                str(self.theta).replace('.', '_'), str(self.mu).replace('.', '_'))
            print("\twriting top %d scores to %s" % (self.num_pred_to_write, out_file))

            with open(out_file, 'w') as out:
                for i in range(Xh.shape[0]):
                    scores = Xh[i].toarray().flatten()
                    # convert the nodes back to their names, and make a dictionary out of the scores
                    scores = {self.prots[j]:s for j, s in enumerate(scores)}
                    self.write_scores_to_file(scores, goid=self.goids[i], file_handle=out,
                            num_pred_to_write=self.num_pred_to_write)
        return Xh

    def run_alg_with_params(self, alg, out_pref=None):
        """ Call the appropriate algorithm's function
        """
        params_results = {}

        # TODO run these separately because they have different parameters
        if alg in ["sinksourceplus", "sinksource", "localplus", "local"]:
            goid_scores, params_results = self.run_ss_with_params(
                    alg, out_pref=out_pref)
        elif alg in ["sinksource-ripple", "sinksourceplus-ripple", 
                'sinksource-squeeze', 'sinksourceplus-squeeze']:
            goid_scores, params_results = self.run_topk_with_params(
                    alg, out_pref=out_pref)
                    #out_pref=out_pref, forced=forced)
        elif alg in ["genemania"]:
            out_file = "%sresults.txt" % (out_pref) if out_pref is not None else None
            goid_scores, params_results = self.run_alg_on_goterms(alg,
                    out_file=out_file)

        return goid_scores, params_results

    def run_ss_with_params(self, alg, out_pref=None):
        all_params_results = {} 
        #for a, eps in tqdm(list(itertools.product(*[a_list, eps_list]))):
        a = self.alpha
        out_file = None
        if out_pref is not None:
            out_file = "%sa%s-eps%s.txt" % (
                    out_pref, str(a).replace('.', '_'), str(self.eps).replace('.', '_'))
            #header = "#%s\tGO term: %s\t%d positives\t%d negatives\ta=%s\teps=%s\n" \
            #        % (alg, goterm, len(positives), len(negatives), str(a), str(eps))

        goid_scores, params_results = self.run_alg_on_goterms(alg,
                out_file=out_file, a=a, eps=self.eps)
        all_params_results.update(params_results)

        if self.verbose:
            print(params_results_to_table(params_results))
        return goid_scores, all_params_results

    def run_topk_with_params(self, alg, out_pref=None):
        all_params_results = {}

        a = self.alpha
        # Run using all combinations of parameters
        # s and t won't change the results, just the amount of time, so loop through those separately
        total_params_list = list(itertools.product(*[self.k_list, self.t_list, self.s_list, self.epsUB_list]))
        #params_list = list(itertools.product(*[k_list, a_list, deltaUBLB_list]))
        #t_s_list = list(itertools.product(*[t_list, s_list]))
        print("Running %d combinations of parameters" % (len(total_params_list)))
        for k, epsUB in itertools.product(*[self.k_list, self.epsUB_list]):

            if "ripple" in alg:
                for t, s in itertools.product(*[self.t_list, self.s_list]):
                    out_file = None
                    if out_pref is not None:
                        out_file = "%sk%d-t%d-s%d-a%s-eps%s-scipy.txt" % (out_pref, k, t, s, 
                                str(a).replace('.', '_'), str(epsUB).replace('.', '_'))

                    #print("Running %s with %d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, deltaUBLB=%s for GO term %s" \
                    #        % (alg, len(positives), len(negatives), k, t, s, str(a), str(deltaUBLB), goterm))

                    goid_scores, params_results = self.run_alg_on_goterms(alg, out_file=out_file,
                            a=a, k=k, t=t, s=s, epsUB=epsUB)
                    all_params_results.update(params_results)
                        # also keep track of the time it takes for each of the parameter sets

            elif 'squeeze' in alg:
                out_file = None
                if out_pref is not None:
                    out_file = "%sk%d-a%s-epsUB%s.txt" % (out_pref, k, str(a).replace('.', '_'),
                                                          str(epsUB).replace('.', '_'))
                goid_scores, params_results = self.run_alg_on_goterms(alg, out_file=out_file,
                        a=a, k=k, epsUB=epsUB)
                all_params_results.update(params_results)

            if self.verbose:
                print(params_results_to_table(all_params_results))

        return goid_scores, all_params_results

    def run_alg_on_goterms(self, alg, out_file=None, 
            a=0.8, eps='-', k='-', t='-', s='-', epsUB=0):
        """ Run the specified algorithm with the given parameters for each goterm 
        *returns*: a dictionary of scores from the algorithm for each goterm
            and a dictionary of summary statistics about the run
        """
        # scores from the algorithm for each goterm
        params_results = {} 
        #goid_scores = {}
        # store the results in a sparse matrix
        goid_scores = sparse.lil_matrix(self.ann_matrix.shape, dtype=np.float)
        try:
            if out_file is not None:
                if self.forcealg is False and os.path.isfile(out_file):
                    print("%s already exists. Use --forcealg to overwrite" % (out_file))
                    print("\tskipping %s" % (alg))
                    return {}, {}

                # saves a lot of time keeping the file handle open
                file_handle = open(out_file, 'w')
                file_handle.write("#goterm\tprot\tscore\n")
                #file_handle.write("#%s\tk=%d\tt=%d\ts=%d\ta=%s\tepsUB=%s\n" \
                #        % (alg, k, t, s, str(a), str(epsUB)))
            if self.compare_ranks and 'squeeze' in alg:
                out_dir = "outputs/viz/ranks"
                utils.checkDir(out_dir)
                k_str = str(k) if self.rank_topk is True else 'all'
                rank_file = "%s/compare-ranks-%s-k%s-%s.txt" % (out_dir, alg, k_str, self.exp_name)
                # write the ranks to a file
                print("Writing rank comparison to: %s" % (rank_file))
                rank_fh = open(rank_file, 'w', buffering=100)
                rank_fh.write("#goterm\tnum_pos\titer\tkendalltau\tspearmanr\tnum_unranked\tmax_d\n")

            print("Running %s for %d goterms. Writing to %s" % (alg, self.ann_matrix.shape[0], out_file))
            for i in tqdm(range(self.ann_matrix.shape[0]), total=self.ann_matrix.shape[0], disable=not self.progress_bar):
                goid = self.goids[i]
                positives, negatives = alg_utils.get_goid_pos_neg(self.ann_matrix, i)
                # get the row corresponding to the current goids annotations 
                if len(positives) < 10:
                    tqdm.write("Skipping goterm %s. It has %d annotations which is < the minimum 10." % (goid, len(positives)))
                    continue
                
                if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[self.version] and not self.unweighted:
                    y = np.zeros(self.normalized_nets[0].shape[0])
                    y[positives] = 1
                    y[negatives] = -1
                    # weight the network for each GO term individually
                    W = setup.weight_net_goterm(y, self.normalized_nets, self.net_names, goid)
                    self.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=self.ss_lambda)
                    if alg == 'genemania':
                        self.L = genemania.setup_laplacian(W)

                scores_arr, curr_params_results, _ = self.run_alg(alg, positives, negatives, 
                        a=a, eps=eps, k=k, t=t, s=s, epsUB=epsUB, goid=goid)
                # storing all of the scores for each goterm takes a lot of memory
                # rather than store the scores in a dictionary, store them in a sparse matrix
                # split the dictionary into a list of indices and a list of scores
                #indices, score_list = zip(*scores.items()) 
                ## build an array of the scores and set it in the goid sparse matrix of scores
                #scores_arr = np.zeros(goid_scores.shape[1])
                #scores_arr[list(indices)] = list(score_list)
                goid_scores[i] = scores_arr
                params_results.update(curr_params_results)

                if self.compare_ranks and 'squeeze' in alg:
                    # compare how long it takes for the ranks to match the previous run
                    tqdm.write("\tRepeating the run, but comparing the ranks from the previous run at each iteration")
                    ranks = [n for n in sorted(scores, key=scores.get, reverse=True)]
                    ranks = ranks[:k] if self.rank_topk is True else ranks
                    _, _, ss_squeeze = self.run_alg(alg, positives, negatives,
                            a=a, eps=eps, k=k, t=t, s=s, epsUB=epsUB, goid=goid, ranks_to_compare=ranks)
                    rank_fh.write(''.join("%s\t%d\t%d\t%0.6f\t%0.6f\t%d\t%s\n" % (goid, len(positives), i+1,
                        ss_squeeze.kendalltau_list[i], ss_squeeze.spearmanr_list[i],
                        ss_squeeze.num_unranked_list[i], ss_squeeze.max_d_list[i])
                                        for i in range(ss_squeeze.num_iters)))

                if out_file is not None:
                    # convert the nodes back to their names, and make a dictionary out of the scores
                    scores = {self.prots[i]:s for i, s in enumerate(scores_arr)}
                    self.write_scores_to_file(scores, goid=goid, file_handle=file_handle,
                            num_pred_to_write=self.num_pred_to_write)

                # TODO move this somewhere else
                # plot the max_ds
                #self.plot_max_ds(max_ds)
        except:
            if out_file is not None:
                file_handle.close()
            raise

        if out_file is not None:
            file_handle.close()
            print("Finished running %s for %d goterms. Wrote to %s" % (alg, len(self.goids), out_file))
        if self.compare_ranks and 'squeeze' in alg:
            rank_fh.close()
            print("Finished running rank comparison. Wrote to: %s" % (rank_file))

        return goid_scores, params_results

    def run_alg(self, alg, positives, negatives, nodes_to_rank=None,
                a=0.8, eps='-', k='-', t='-', s='-', epsUB=0, goid='-', ranks_to_compare=None):
        """ Run the specified algorithm with the given parameters for each goterm 
        *returns*: a dictionary of scores from the algorithm for each goterm
            and a dictionary of summary statistics about the run
        """
        num_unk = self.P.shape[0] - len(positives)
        if alg in POSITIVE_ALGS:
            negatives=None
        else:
            num_unk -= len(negatives) 
        len_N = self.P.shape[0]
        ss_obj = None 

        # TODO streamline calling the correct function. They all take the same parameters
        # This uses the same UB as Ripple
        if alg in ['sinksourceplus', 'sinksource']:
            scores, total_time, iters, comp = sinksource.runSinkSource(
                    self.P, positives, negatives=negatives, max_iters=self.max_iters, delta=eps, a=a)
            update_time = '-'
        elif alg in ['localplus', 'local']:
            scores, total_time = sinksource.runLocal(
                    self.P, positives, negatives=negatives)
            update_time, total_time, iters, comp = [-1, total_time, 1, -1]
        elif alg == 'genemania':
            y = np.zeros(self.P.shape[0])
            y[positives] = 1
            y[negatives] = -1
            scores, total_time = genemania.runGeneMANIA(self.L, y)
            update_time, iters, comp = [-1, -1, -1]
        elif alg in ['sinksourceplus-squeeze', 'sinksource-squeeze']:
            ss_obj = sinksource_squeeze.SinkSourceSqueeze(
                    self.P, positives, negatives=negatives, k=k, a=a, epsUB=epsUB, verbose=self.verbose,
                    rank_topk=self.rank_topk, rank_all=self.rank_all, rank_nodes=nodes_to_rank,
                    ranks_to_compare=ranks_to_compare)
            R, scores = ss_obj.runSinkSourceSqueeze() 
            total_time, update_time, iters, comp, len_N = ss_obj.get_stats()
            num_unk = len_N
        elif alg in ['sinksourceplus-ripple', 'sinksource-ripple']:
            ss_obj = sinksource_ripple.SinkSourceRipple(
                    self.P, positives, negatives=negatives, k=k, a=a, t=t, s=s, epsUB=epsUB,
                    verbose=self.verbose)
            R, scores = ss_obj.runSinkSourceRipple() 
            total_time, update_time, iters, comp, len_N = ss_obj.get_stats()
            num_unk = ss_obj.P.shape[0]

        tqdm.write("\t%s converged after %d iterations " % (alg, iters) +
                "(%0.2f sec) for goterm %s" % (total_time, goid))

        # also keep track of the time it takes for each of the parameter sets
        params_results = {} 
        params_key = (alg, goid, len(positives), num_unk, a, eps, k, t, s, epsUB)
        params_results[params_key] = (total_time, update_time, iters, comp, len_N)

        # TODO return a reference to the object used to run the algorithm to get additional statistics where necessary
        return scores, params_results, ss_obj

    def plot_max_ds(self, max_d_list):
#import matplotlib.pyplot as plt
#plt.plot(max_d_list)
#plt.yaxis("Maximum score change")
#plt.xaxis("iteration")
        out_file = "outputs/viz/max-d.txt"
        print("writing %s" % out_file)
        with open(out_file, 'w') as out:
            out.write('\n'.join(str(x) for x in max_d_list))
        #plt.savefig(out_file)

    def write_scores_to_file(self, scores, goid='', out_file=None, file_handle=None,
            num_pred_to_write=100, header="", append=True):
        """
        *scores*: dictionary of node_name: score
        *num_pred_to_write*: number of predictions (node scores) to write to the file 
            (sorted by score in decreasing order). If -1, all will be written
        """

        if num_pred_to_write == -1:
            num_pred_to_write = len(scores) 

        if out_file is not None:
            if append:
                print("Appending %d scores for goterm %s to %s" % (num_pred_to_write, goid, out_file))
                out_type = 'a'
            else:
                print("Writing %d scores for goterm %s to %s" % (num_pred_to_write, goid, out_file))
                out_type = 'w'

            file_handle = open(out_file, out_type)
        elif file_handle is None:
            print("Warning: out_file and file_handle are None. Not writing scores to a file")
            return 

        # write the scores to a file
        # for now, write the top k node's results
        # TODO Possibly write the results to a JSON file so it's easier to load
        file_handle.write(header)
        for n in sorted(scores, key=scores.get, reverse=True)[:num_pred_to_write]:
            file_handle.write("%s\t%s\t%s\n" % (goid, n, str(scores[n])))
        return

    def run_goterm_cv(self, alg, folds=5, out_pref=None):
        """
        For every GO term, split the positives and negatives into the specified number of folds 
            and run and evaluate the given algorithm
        """
        # compare how long it takes to run each from scratch vs using the previous run's scores
        params_results = {}
        #overall_time_normal = 0
        #overall_time_using_predictions = 0

        #if alg not in ['sinksource', "sinksourceplus", "sinksource-squeeze", "sinksourceplus-squeeze"]:
        #    print("%s not yet implemented for CV" % (alg))
        #    return params_results
    #    if alg == "sinksourceplus":
    #        # this should already have been created from the regular run
    #        out_file = "%sa%s-eps%s.txt" % (out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'))
    #        print("\tReading scores from %s" % (out_file))
    #        prediction_scores = {node2idx[goterm]: float(score) for goterm, score in utils.readColumns(out_file, 1, 2)}
    #        #predition_scores, time, iters = sinksource.runSinkSource(H, positives, negatives=None, max_iters=1000, delta=eps, a=a)

        print("Running cross-validation for %d goterms using %d folds for %s; a=%s, eps=%s" % (
            self.ann_matrix.shape[0], self.cross_validation_folds, alg, self.alpha, self.eps))

        out_file = None
        if out_pref is not None:
            out_file = "%sa%s-eps%s.txt" % (
                    out_pref, str(self.alpha).replace('.', '_'), str(self.eps).replace('.', '_'))
            print("Writing CV results to %s" % (out_file))
            file_handle = open(out_file, 'w')
            file_handle.write("#goterm\tfmax\n")

        #for goterm in tqdm(goid_pos_neg):
        #    positives, negatives = goid_pos_neg[goterm]['pos'], goid_pos_neg[goterm]['neg']
        for i in tqdm(range(self.ann_matrix.shape[0]), total=self.ann_matrix.shape[0], disable=not self.progress_bar):
            goid = self.goids[i]
            positives, negatives = alg_utils.get_goid_pos_neg(self.ann_matrix, i)
            print("%d positives, %d negatives for goterm %s" % (len(positives), len(negatives), goid))
            kf = KFold(n_splits=folds, shuffle=True)
            kf_neg = KFold(n_splits=folds, shuffle=True)
            kf.get_n_splits(positives)
            kf_neg.get_n_splits(negatives)
            fold = 0
            # because each fold contains a different set of positives, and combined they contain all positives,
            # store all of the prediction scores from each fold in a list
            combined_fold_scores = {}
            if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[self.version] and not self.unweighted:
                unknowns = set(range(self.normalized_nets[0].shape[0])) - set(positives) 
            else:
                unknowns = set(range(self.P.shape[0])) - set(positives) 
            if 'plus' not in alg:
                unknowns = unknowns - set(negatives)

            for (pos_train_idx, pos_test_idx), (neg_train_idx, neg_test_idx) in zip(kf.split(positives), kf_neg.split(negatives)):
                fold += 1
                pos_train, pos_test = positives[pos_train_idx], positives[pos_test_idx]
                neg_train, neg_test = negatives[neg_train_idx], negatives[neg_test_idx]

                #if self.rank_pos is True and ('squeeze' in alg or 'ripple' in alg):
                nodes_to_rank = None 
                if self.rank_pos is True:
                    nodes_to_rank = set(list(pos_test))

                # 
                if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[self.version] and not self.unweighted:
                    y = np.zeros(self.normalized_nets[0].shape[0])
                    y[pos_train] = 1
                    y[neg_train] = -1

                    #W = setup.weight_SWSN(curr_train_ann_mat, self.sparse_networks,
                    #        net_names=self.net_names, out_file=out_file, nodes=prots)
                    # TODO add an option for this
                    # for now, weight the network for each GO term individually
                    W = setup.weight_net_goterm(y, self.normalized_nets, self.net_names, goid)
                    self.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=self.ss_lambda)
                    if alg == 'genemania':
                        self.L = genemania.setup_laplacian(W)

                # TODO use the list of parameters (e.g., alpha, eps) rather than the first one.
                scores, params_results, _ = self.run_alg(alg, pos_train, neg_train,
                        nodes_to_rank=nodes_to_rank, goid=goid,
                        a=self.alpha, eps=self.eps, k=self.k_list[0],
                        t=self.t_list[0], s=self.s_list[0], epsUB=self.epsUB_list[0])
                    # scores is a dictionary of node integers 
                    # containing only scores for the non-positive and non-negative nodes
                    #scores = np.array([scores[n] for n in pos_test])
                    #test_scores[test_idx] = scores 

                # add a check to ensure the scores are actually available
                if len(scores) == 0:
                    print("WARNING: No scores found. Skipping")
                    continue

                # the test positives and negatives will appear in a single fold
                fold_scores = {n:scores[n] for n in set(pos_test).union(set(neg_test))}
                if nodes_to_rank is not None:
                    fold_pos_ranks = [i for i, x in enumerate(sorted(scores, key=scores.get)) if x in nodes_to_rank]
                    print("Ranks of left out positives:")
                    print(fold_pos_ranks)
                # the unknowns will be in each fold, so append the fold number to those nodes
                #for n in unknowns:
                #    fold_scores["%d-%d" % (n, fold)] = scores[n]
                combined_fold_scores.update(fold_scores)

                #prec, tpr, fpr = self.compute_eval_measures(fold_scores, pos_test, neg_test)
                #fmax = self.compute_fmax(prec, tpr)
                #print("\tfold %d fmax: %0.3f" % (fold, fmax))

            if len(combined_fold_scores) == 0:
                continue

            # sort the combined scores by the score, and then compute the metrics on the combined scores
            prec, recall, fpr = self.compute_eval_measures(combined_fold_scores, positives, negatives)
            fmax = self.compute_fmax(prec, recall)
            tqdm.write("\toverall fmax: %0.3f" % (fmax))
            if out_file is not None:
                file_handle.write("%s\t%0.4f\n" % (goid, fmax))
#            params_key = (alg, goterm, 'folds=5', '-', '-', a, 'eps=%s'%str(eps))
#            params_results[params_key] = (time_diff, avg_iter_diff, '-')
#            params_key = (alg, goterm, len(positives), num_unk, a, eps, k, t, s, epsUB)
#            params_results[params_key] = (time, iters, comp, len_N)
        if out_file is not None:
            file_handle.close()
            print("Finished running %s for %d goterms. Wrote to %s" % (alg, self.ann_matrix.shape[0], out_file))

        return params_results

    def evaluate_ground_truth(
            self, goid_scores, goids, true_ann_matrix, out_pref,
            non_pos_as_neg_eval=False, taxon='-',
            alg='', write_prec_rec=False):

        score_goids2idx = {g: i for i, g in enumerate(goids)}
        print("Computing fmax from ground truth of %d goterms" % (true_ann_matrix.shape[0]))
        goid_stats = {}
        goid_num_pos = {} 
        goid_prec_rec = {}
        #curr_goid2idx = {g: i for i, g in enumerate(goids)}
        for i in range(true_ann_matrix.shape[0]):
            goid = self.goids[i]
            # make sure the scores are actually available first
            #if goid not in goid_scores:
            #    print("WARNING: goid %s not in goid_scores" % (goid))
            #    continue
            if goid not in self.goid2idx:
                print("WARNING: goid %s not in initial set of %d goids" % (
                                  goid, len(self.goids)))
                continue
            #positives, negatives = goid_pos_neg[goterm]['pos'], goid_pos_neg[goterm]['neg']
            # get the row corresponding to the current goids annotations 
            goid_ann = true_ann_matrix[i,:].toarray().flatten()
            positives = np.where(goid_ann > 0)[0]
            # to get the scores, map the current goid index to the
            # index of the goid in the scores matrix
            scores = goid_scores[score_goids2idx[goid]].toarray().flatten()
            # convert the scores to a dictionary
            scores = {i: s for i,s in enumerate(scores)}
            # this was already done
            #positives = set(self.node2idx[n] for n in positives if n in self.node2idx)
            goid_num_pos[goid] = len(positives)
            if len(positives) == 0:
                print("%s has 0 positives after restricting to nodes in the network. Skipping" % (goid))
                continue
            if non_pos_as_neg_eval is True:
                # leave everything not a positive as a negative
                negatives = None
            else:
                # alternatively, use the negatives from that species as the set of negatives
                negatives = np.where(goid_ann < 0)[0]
                if len(negatives) == 0:
                    print("WARNING: 0 negatives for %s - %s. Skipping" % (goid, taxon))
                    continue
            prec, recall, fpr, pos_neg_stats = self.compute_eval_measures(scores, positives, negatives=negatives, track_pos_neg=True)
            # TODO only write the prec and rec and the points that change in recall(?)
            # TODO also write the prec, recall
            if write_prec_rec:
                goid_prec_rec[goid] = (prec, recall, pos_neg_stats)
            #print((len(scores), len(positives), len(prec), len(recall)))
            fmax = self.compute_fmax(prec, recall)
            avgp = self.compute_avgp(prec, recall)
            if len(prec) == 1:
                auprc = 0
                auroc = 0
            else:
                auprc = self.compute_auprc(prec, recall)
                auroc = self.compute_auroc([r for r, f in fpr], [f for r, f in fpr])
            goid_stats[goid] = (fmax, avgp, auprc, auroc)
            if self.verbose:
                print("%s fmax: %0.4f" % (goid, fmax))

        if not write_prec_rec:
            # TODO
            #out_file = "%s.txt" % (out_pref)
            #if alg in ['sinksource', 'sinksourceplus']:
            out_file = "%sa%s.txt" % (
                    out_pref, str(self.alpha).replace('.', '_'))
            # don't write the header each time
            if not os.path.isfile(out_file):
                print("Writing results to %s" % (out_file))
                with open(out_file, 'w') as out:
                    if taxon == '-':
                        out.write("#goid\tfmax\tavgp\tauprc\tauroc\t# ann\n")
                    else:
                        out.write("#taxon\tgoid\tfmax\tavgp\tauprc\tauroc\t# ann\n")
            else:
                print("Appending results to %s" % (out_file))
            with open(out_file, 'a') as out:
                out.write(''.join(["%s%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%d\n" % (
                    "%s\t"%taxon if taxon != "-" else "",
                    g, fmax, avgp, auprc, auroc, goid_num_pos[g]
                    ) for g, (fmax, avgp, auprc, auroc) in goid_stats.items()]))

        if write_prec_rec:
            out_file = "%sprec-rec%s%s.txt" % (
                out_pref, taxon, '-%s'%(list(goid_prec_rec.keys())[0]) if len(goid_prec_rec) == 1 else "")
            if alg in ['sinksource', 'sinksourceplus']:
                out_file = "%sa%s-prec-rec%s%s.txt" % (
                    out_pref, str(self.alpha).replace('.', '_'), taxon,
                    '-%s'%(list(goid_prec_rec.keys())[0]) if len(goid_prec_rec) == 1 else "")
            print("writing prec/rec to %s" % (out_file))
            with open(out_file, 'w') as out:
                out.write("goid\tprec\trec\tnode\tscore\tidxpos/neg\n")
                for goid, (prec, rec, pos_neg_stats) in goid_prec_rec.items():
                    out.write(''.join(["%s\t%0.4f\t%0.4f\t%s\t%0.4f\t%d\t%d\n" % (
                        goid, p, r, self.prots[n], s, idx, pos_neg) for p,r,(n,s,idx,pos_neg) in zip(prec[1:], rec[1:], pos_neg_stats)]))

    def compute_eval_measures(self, scores, positives, negatives=None, track_pos_neg=False):
        """
        Compute the precision and false-positive rate at each change in recall (true-positive rate)
        *scores*: dictionary containing a score for each node
        *negatives*: if negatives are given, then the FP will only be from the set of negatives given
        *track_pos_neg*: if specified, track the score and rank of the positive and negative nodes,
            and return a tuple of the node ids in order of their score, their score, their idx, and 1/-1 for pos/neg
        """
        #f1_score = metrics.f1score(positives, 
        #num_unknowns = len(scores) - len(positives) 
        positives = set(positives)
        check_negatives = False
        if negatives is not None:
            check_negatives = True 
            negatives = set(negatives)
        else:
            print("TODO. Treating all non-positives as negatives not yet implemented.")
        # compute the precision and recall at each change in recall
        nodes_sorted_by_scores = sorted(scores, key=scores.get, reverse=True)
        #print("computing the rank of positive nodes")
        # this is really slow...
        #pos_ranks = sorted([nodes_sorted_by_scores.index(p)+1 for p in positives])
        #print("%d positives, %d pos_ranks" % (len(positives), len(pos_ranks)))
        #print(pos_ranks)
        #print([scores[s] for s in nodes_sorted_by_scores[:pos_ranks[0]+1]])
        precision = [1]
        recall = [0]
        fpr = []
        pos_neg_stats = []  # tuple containing the node, score and idx
        # TP is the # of correctly predicted positives so far
        TP = 0
        FP = 0
        rec = 0
        for i, n in enumerate(nodes_sorted_by_scores):
            # TODO this could be slow if there are many positives
            if n in positives:
                TP += 1
                # precisions is the # of true positives / # true positives + # of false positives (or the total # of predictions)
                precision.append(TP / float(TP + FP))
                # recall is the # of recovered positives TP / TP + FN (total # of positives)
                rec = TP / float(len(positives))
                recall.append(rec)
                # fpr is the FP / FP + TN
                fpr.append((rec, FP / float(len(negatives))))
                if track_pos_neg:
                    pos_neg_stats.append((n, scores[n], i, 1)) 
            elif check_negatives is False or n in negatives:
                FP += 1
                fpr.append((rec, FP / float(len(negatives))))
            #else:
            #    continue

        # TODO how should I handle this case?
        if len(precision) == 0:
            precision.append(0)
            recall.append(1)

        #print(precision[0], recall[0], fpr[0])

        if track_pos_neg:
            return precision, recall, fpr, pos_neg_stats
        else:
            return precision, recall, fpr
     
    def compute_fmax(self, prec, rec):
        f_measures = []
        for i in range(len(prec)):
            p, r = prec[i], rec[i]
            if p+r == 0:
                harmonic_mean = 0
            else:
                harmonic_mean = (2*p*r)/(p+r)
            f_measures.append(harmonic_mean)
        return max(f_measures)

    def compute_avgp(self, prec, rec):
        # average precision score
        # see http://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html#sklearn.metrics.average_precision_score
        avgp = 0
        prev_r = 0 
        for p,r in zip(prec, rec):
            recall_change = r - prev_r
            avgp += (recall_change*p)
            prev_r = r
        #avgp = avgp / float(len(alg_prec_rec))
        return avgp

    def compute_auprc(self, prec, rec):
        auprc = metrics.auc(rec, prec)
        return auprc

    def compute_auroc(self, tpr, fpr):
        auroc = metrics.auc(fpr, tpr)
        return auroc


def params_results_to_table(params_results):
    # print out a table of the results for each parameter set
    # print out the table after each iteration so the results are still viewable if the script quits earl or something
    results = "\t".join(['alg', 'goterm', '# pos', '# unk', 'a', 'k', 't', 's', 'eps', 'epsUB', 'time', 'update-time', 'iters', '# comp', 'len_N']) + '\n'
    for params_key in sorted(params_results):
        alg, goterm, num_pos, num_unk, a, eps, k, t, s, epsUB = params_key
        time, update_time, iters, total_comp, len_N = params_results[params_key]
        update_time = "%0.3f" % update_time if isinstance(update_time, float) else update_time
        results += "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%0.3f\t%s\t%d\t%0.2e\t%s\n" % \
                (alg, goterm, num_pos, num_unk, str(a), str(k), str(t), str(s), str(eps), str(epsUB), 
                        time, update_time, iters, total_comp, str(len_N))
    return results


def parse_args(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    # general parameters
    group = OptionGroup(parser, 'Main Options')
    group.add_option('','--version', type='string', default="2017_10-seq-sim",
                     help="Version of the PPI to run. Default: %s" % ("2017_10-seq-sim") + \
                     "\nOptions are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    group.add_option('-N','--net-file', type='string',
                     help="Network file to use. Default is the version's default network")
    group.add_option('-A', '--algorithm', action="append",
                     help="Algorithm for which to get predictions. Default is all of them. Options: '%s'" % ("', '".join(ALGORITHMS)))
    group.add_option('', '--exp-name', type='string',
                     help="Experiment name to use when running GAIN.")
    group.add_option('', '--pos-neg-file', type='string', action='append',
                     help="File containing positive and negative examples for each GO term")
    group.add_option('', '--only-functions', type='string',
                     help="Run GAIN using only the functions in a specified file (should be in XX format i.e., without the leading GO:00).")
    group.add_option('-G', '--goterm', type='string', action="append",
                     help="Specify the GO terms to use (should be in GO:00XX format)")
    parser.add_option_group(group)

    # parameters for running algorithms
    group = OptionGroup(parser, 'Algorithm options')
    group.add_option('', '--weight-swsn', action="store_true", default=False,
                     help="Option to integrate multiple networks using the SWSN (Simultaneous Weighting with Specific Negatives) method. Only takes effect if multiple networks are present for the specified version. Default=False (integrate the networks for each GO term individually).")
    group.add_option('', '--unweighted', action="store_true", default=False,
                     help="Option to ignore edge weights when running algorithms. Default=False (weighted)")
    group.add_option('-l', '--sinksourceplus-lambda', type=float, 
                     help="lambda parameter to specify the weight connecting the unknowns to the negative 'ground' node. Default=None")
    group.add_option('-a', '--alpha', type=float, default=0.8,
                     help="Alpha insulation parameter. Default=0.8")
    group.add_option('', '--eps', type=float, default=0.0001,
                     help="Stopping criteria for SinkSource. Default=0.0001")
    group.add_option('', '--max-iters', type=int, default=1000,
                     help="Maximum # of iterations for SinkSource. Default=1000")
    # ripple/squeeze parameters
    group.add_option('-k', '--k', type=int, action="append",
                     help="Top-k for Ripple, Default=200")
    group.add_option('-e', '--epsUB', type=float, action="append",
                     help="Parameter to return the top-k if all other nodes have an UB - epsUB < the kth node's LB. Default=0")
    group.add_option('-t', '--t', type=int, action="append",
                     help="t parameter for Ripple. Default=2")
    group.add_option('-s', '--s', type=int, action="append",
                     help="s parameter for Ripple. Default=200")
    group.add_option('', '--rank-topk', action="store_true", default=False,
                     help="Continue iterating until the top-k nodes ranks are fixed (comparing the UB and LB). Currently only available for SinkSourceSqueeze")
    group.add_option('', '--rank-all', action="store_true", default=False,
                     help="Continue iterating until all nodes ranks are fixed (comparing the UB and LB). Currently only available for SinkSourceSqueeze")
    group.add_option('', '--rank-pos', action="store_true", default=False,
                     help="During cross-validation, continue iterating until the nodes of the left out positives are fixed (comparing the UB and LB). Currently only available for SinkSourceSqueeze")
    group.add_option('', '--compare-ranks', action="store_true", default=False,
                     help="Compare how long it takes (# iters) for the ranks to match the final fixed ranks." +
                     "Currently only implemented with ss_squeeze and --rank-topk, --rank-all")
    parser.add_option_group(group)

    # birgrank options
    group = OptionGroup(parser, 'BirgRank options')
    group.add_option('-b', '--obo-file', type='string', default=f_settings.GO_FILE,
                     help="GO OBO file which contains the GO DAG. Used if running AptRank/BirgRank. Default: %s" % (f_settings.GO_FILE))
    group.add_option('', '--theta', type=float, default=.5,
                     help="BirgRank parameter: (1-theta) percent of Rtrain used in seeding vectors")
    group.add_option('', '--mu', type=float, default=.5,
                     help="BirgRank parameter: (1-mu) percent of random walkers diffuse from G via Rtrain to H")
    parser.add_option_group(group)

    # parameters for STRING networks
    group = OptionGroup(parser, 'STRING options')
    group.add_option('', '--string-combined', action="store_true", default=False,
            help="Use only the STRING combined network: \n\tcombined_score")
    group.add_option('', '--string-core', action="store_true", default=False,
            help="Use only the 6 core networks: \n\t%s" % (', '.join(setup.CORE_STRING_NETWORKS)))
    group.add_option('', '--string-non-transferred', action="store_true", default=False,
            help="Use all non-transferred networks: \n\t%s" % (', '.join(setup.NON_TRANSFERRED_STRING_NETWORKS)))
    group.add_option('', '--string-all', action="store_true", default=False,
            help="Use all individual 13 STRING networks: \n\t%s" % (', '.join(setup.STRING_NETWORKS)))
    group.add_option('-S', '--string-networks', type='string', 
            help="Comma-separated list of string networks to use. " +
                 "If specified, other STRING options will be ignored." +
                 "Default: %s" % (', '.join(setup.CORE_STRING_NETWORKS)))
    parser.add_option_group(group)

    # additional parameters
    group = OptionGroup(parser, 'Additional options')
    group.add_option('-W', '--num-pred-to-write', type='int', default=100,
                     help="Number of predictions to write to the file. If 0, none will be written. If -1, all will be written. Default=100")
    group.add_option('', '--only-cv', action="store_true", default=False,
                     help="Perform cross-validation only")
    group.add_option('-C', '--cross-validation-folds', type='int',
                     help="Perform cross validation using the specified # of folds. Usually 5")
    # TODO finish adding this option
    #group.add_option('-T', '--ground-truth-file', type='string',
    #                 help="File containing true annotations with which to evaluate predictions")
    group.add_option('', '--forcealg', action="store_true", default=False,
                     help="Force re-running algorithms if the output files already exist")
    group.add_option('', '--forcenet', action="store_true", default=False,
                     help="Force re-building network matrix from scratch")
    group.add_option('', '--verbose', action="store_true", default=False,
                     help="Print additional info about running times and such")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args(args)

    if opts.exp_name is None or opts.pos_neg_file is None:
        print("--exp-name, --pos-neg-file, required")
        sys.exit(1)

    # if neither are provided, just use all GO terms in the pos/neg file
#    if opts.goterm is None and opts.only_functions is None:
#        print("--goterm or --only_functions required")
#        sys.exit(1)

    if opts.algorithm is None:
        opts.algorithm = ALGORITHMS

    if opts.version not in f_settings.ALLOWEDVERSIONS:
        print("ERROR: '%s' not an allowed version. Options are: %s." % (opts.version, ', '.join(f_settings.ALLOWEDVERSIONS)))
        sys.exit(1)

    #if opts.algorithm not in f_settings.ALGORITHM_OPTIONS:
    #    print "--algorithm %s is not a valid algorithm name" % (opts.algorithm)
    #    sys.exit(1)

    # TODO
    for alg in opts.algorithm:
        if alg not in ALGORITHMS:
            print("ERROR: '%s' not a valid algorithm name. Algorithms are: '%s'." % (alg, ', '.join(ALGORITHMS)))
            sys.exit(1)

    opts.k = opts.k if opts.k is not None else [200]
    opts.t = opts.t if opts.t is not None else [2]
    opts.s = opts.s if opts.s is not None else [200]
    opts.epsUB = opts.epsUB if opts.epsUB is not None else [0]

    # setup the selection of string networks 
    string_networks = []
    if opts.string_networks:
        string_networks = opts.string_networks.split(',')
        for net in string_networks:
            if net not in setup.STRING_NETWORKS:
                print("ERROR: STRING network '%s' not one of the" +
                      "available choices which are: \n\t%s" % (net, ', '.join(setup.STRING_NETWORKS)))
                sys.exit(1)
    elif opts.string_combined:
        string_networks = ['combined_score']
    elif opts.string_core:
        string_networks = setup.CORE_STRING_NETWORKS
    elif opts.string_non_transferred:
        string_networks = setup.NON_TRANSFERRED_STRING_NETWORKS
    elif opts.string_all:
        string_networks = setup.STRING_NETWORKS
    else:
        string_networks = setup.CORE_STRING_NETWORKS
    opts.string_networks = string_networks

    return opts


def run():
    #versions = ["2017_10-seq-sim", "2017_10-seq-sim-x5-string"]
    opts = parse_args(sys.argv)
    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 

    #goid_pos, goid_neg = parse_pos_neg_files(opts.pos_neg_file) 
    # load the network matrix and protein IDs
    net_file = opts.net_file
    if net_file is None:
    #    _, _, net_file, _ = f_settings.set_version(opts.version) 
    #else:
        INPUTSPREFIX, _, net_file, selected_strains = f_settings.set_version(opts.version) 
    # TODO this should be better organized so that any STRING networks
    # can be used
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] and not opts.unweighted:
        out_pref_net = "%s/sparse-nets/" % (INPUTSPREFIX)
        utils.checkDir(out_pref_net)
        # build the file containing the sparse networks
        sparse_networks, network_names, prots = setup.create_sparse_net_file(
            opts.version, out_pref_net, selected_strains=selected_strains,
            string_nets=opts.string_networks, string_cutoff=f_settings.STRING_CUTOFF,
            forcenet=False)
        # TODO organize this better
        W = (sparse_networks, network_names)
    else:
        W, prots = alg_utils.setup_sparse_network(net_file)
    # now build the annotation matrix
    ann_matrix, goids = setup.setup_sparse_annotations(opts.pos_neg_file, goterms, prots,
                             selected_species=None, taxon=None)

    # TODO make this more streamlined
    aptrank_data = None 
    if 'birgrank' in opts.algorithm:
        dag_matrix, pos_matrix, dag_goids = run_birgrank.setup_h_ann_matrices(
                prots, opts.obo_file, opts.pos_neg_file, goterms=goids)
        aptrank_data = (dag_matrix, pos_matrix, dag_goids)

    alg_runner = Alg_Runner(
        opts.version, opts.exp_name,
        W, prots, ann_matrix, goids,
        #goid_pos, goid_neg, goterms, opts.net_file, opts.algorithm,
        algorithms=opts.algorithm, weight_swsn=opts.weight_swsn, 
        unweighted=opts.unweighted, ss_lambda=opts.sinksourceplus_lambda,
        alpha=opts.alpha, eps=opts.eps, max_iters=opts.max_iters,
        k_list=opts.k, t_list=opts.t, s_list=opts.s, epsUB_list=opts.epsUB, 
        rank_topk=opts.rank_topk, rank_all=opts.rank_all, rank_pos=opts.rank_pos, compare_ranks=opts.compare_ranks,
        num_pred_to_write=opts.num_pred_to_write, aptrank_data=aptrank_data,
        only_cv=opts.only_cv, cross_validation_folds=opts.cross_validation_folds,
        forcealg=opts.forcealg, forcenet=opts.forcenet, verbose=opts.verbose)
    alg_runner.main()

#    ground_truth_pos = None 
#    if opts.ground_truth_file is not None:
#        print("Ground truth evaluation not yet implemented. Quitting")
#        sys.exit()
#        # TODO implement incorporating negatives
#        ground_truth_pos, ground_truth_neg = alg_utils.parse_pos_neg_file(opts.ground_truth_file)
#        #for goid, prot in utils.readColumns(self.ground_truth_file, 1, 2):
#
##    # TODO also setup the ground truth
##    if opts.ground_truth_file is not None:
##        out_pref = "%s/ground-truth-%sl%d-" % (out_dir, 
##                'unw-' if self.unweighted else '', 0 if self.ss_lambda is None else int(self.ss_lambda))
##        self.evaluate_ground_truth(goid_scores, true_ann_matrix, goids, out_pref,
##                                    non_pos_as_neg_eval=False, taxon=ground_truth_taxon,
##                                    write_prec_rec=write_prec_rec)


if __name__ == "__main__":
    run()
    #main()
