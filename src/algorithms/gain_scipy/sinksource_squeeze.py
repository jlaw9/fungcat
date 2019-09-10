
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import sys
import time
import operator
from collections import defaultdict
# expects python 3 and networkx 2
#import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
import alg_utils
import random
from scipy.stats import kendalltau  #, spearmanr, weightedtau


class SinkSourceSqueeze:

    def __init__(self, P, positives, negatives=None, k=100, a=0.8, 
                 epsUB=0, rank_topk=False, rank_all=False, rank_pos_neg=None,
                 verbose=False, ranks_to_compare=None, scores_to_compare=None,
                 max_iters=1000):
        """
        *P*: Row-normalized sparse-matrix representation of the graph
        *f*: initial vector f of amount of score received from positive nodes
        *epsUB*: if all other nodes have an UB - epsUB < the kth node's LB, then return the top-k 
        *rank_topk*: require that the ranks of the top-k nodes be fixed using their UB and LB
        *rank_all*: require that the ranks of all nodes be fixed using their UB and LB
        *rank_pos_neg*: tuple of sets of positive and negative nodes. 
            We only require that the LB and UB of positives not overlap with any negative nodes.
            The k parameter will be ignored.
        *ranks_to_compare*: A list of nodes where the index of the node in the list is the rank of that node.
            For example, if node 20 was ranked first and node 50 was ranked second, the list would have [20, 50]
            Used to compare with the current ranking after each iteration.
        *scores_to_compare*: A vector of the scores of each node to compare with the final scores
        *max_iters*: Maximum number of iterations to run power iteration
        *returns*: The set of top-k nodes, and current scores for all nodes
        """
        self.P = P
        self.positives = positives
        self.negatives = negatives
        self.k = k
        self.a = a
        self.epsUB = epsUB
        self.rank_topk = rank_topk
        self.rank_all = rank_all
        self.num_subsample = 500  # num nodes to sample uniformly at random to check if their ranking is fixed
        self.rank_pos_neg = rank_pos_neg
        self.ranks_to_compare = ranks_to_compare
        self.scores_to_compare = scores_to_compare
        self.max_iters = max_iters
        self.verbose = verbose

        # make sure there are not conflicting options specified
        if self.rank_topk and self.rank_all:
            print("ERROR: cannot specify both rank_topk and rank_all. Please choose one. Quitting")
            sys.exit()
        elif self.rank_topk and self.rank_pos_neg is not None:
            print("ERROR: cannot specify both rank_topk and rank_pos_neg. Please choose one. Quitting")
            sys.exit()
        elif self.rank_all and self.rank_pos_neg is not None:
            print("ERROR: cannot specify both rank_all and rank_pos_neg. Please choose one. Quitting")
            sys.exit()

    def runSinkSourceSqueeze(self):
        self.num_nodes = self.P.shape[0]
        # TODO this should be done once before all predictions are being made
        # check to make sure the graph is normalized because making a copy can take a long time
        #G = alg_utils.normalizeGraphEdgeWeights(G)
        self.P, self.f, self.node2idx, self.idx2node = alg_utils.setupScores(
            self.P, self.positives, self.negatives, a=self.a, 
            remove_nonreachable=True, verbose=self.verbose)
        if len(self.f) == 0:
            print("WARNING: no unknown nodes were reachable from a positive (P matrix and f vector empty after removing nonreachable nodes).")
            print("Setting all scores to 0")
            return [], defaultdict(int)
        # if rank_nodes is specified, map those node ids to the current indices
        #if self.rank_nodes is not None:
        #    self.rank_nodes = set(self.node2idx[n] for n in self.rank_nodes if n in self.node2idx)
        if self.rank_pos_neg is not None:
            unr_pos_nodes, unr_neg_nodes = self.rank_pos_neg
            # some of them could've been unreachable, so remove those and fix the mapping
            self.unr_pos_nodes = set(self.node2idx[n] \
                    for n in unr_pos_nodes if n in self.node2idx)
            self.unr_neg_nodes = set(self.node2idx[n] \
                    for n in unr_neg_nodes if n in self.node2idx)

        if self.verbose:
            if self.negatives is not None:
                print("\t%d positives, %d negatives, %d unknowns, k=%d, a=%s, epsUB=%s"
                        % (len(self.positives), len(self.negatives), self.P.shape[0],
                           self.k, str(self.a), str(self.epsUB)))
            else:
                print("\t%d positives, %d unknowns, k=%d, a=%s, epsUB=%s"
                        % (len(self.positives), self.P.shape[0], self.k, str(self.a), str(self.epsUB)))

        #R, all_LBs, overall_time, iters, comp, max_ds = SinkSourceSqueeze(
        R, all_LBs = self._SinkSourceSqueeze()

        # convert the indexes/node integers back to their original node IDs.
        R = set([self.idx2node[n] for n in R])
        # set the default score to 0 for the rank_nodes as some of them may be unreachable from positives
        #scores = defaultdict(int)
        #for n in range(len(all_LBs)):
        #    scores[self.idx2node[n]] = all_LBs[n]
        #scores = {self.idx2node[n]:all_LBs[n] for n in range(len(all_LBs))}
        scores_arr = np.zeros(self.num_nodes)
        indices = [self.idx2node[n] for n in range(len(all_LBs))]
        scores_arr[indices] = all_LBs
        self.scores_to_compare = all_LBs

        if self.verbose:
            print("SinkSourceSqueeze found top k after %d iterations (%0.3f total sec, %0.3f sec to update)"
                  % (self.num_iters, self.total_time, self.total_update_time))

        return R, scores_arr

    def _SinkSourceSqueeze(self):
        """
        *returns*: The set of top-k nodes, and current scores for all nodes
        """
        # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
        R, unranked_nodes, LBs, prev_LBs, UBs = self.initialize_sets()
        # use this variable to indicate if we are ranking a subsample of the nodes first
        initial_unranked_nodes = None
        # comment out the section below to not subsample the unranked nodes
#        if self.rank_all is True:
#            if self.verbose:
#                print("\tRunning until all nodes are ranked correctly")
#            if len(unranked_nodes) > self.num_subsample:
#                # see if sampling a subset of nodes will speed-up the ranking
#                initial_unranked_nodes = unranked_nodes.copy()
#                # randomly sample 1/100 of the nodes to check for fixed ranks
#                #unranked_nodes = random.sample(initial_unranked_nodes, len(unranked_nodes) / float(100))
#                unranked_nodes = set(random.sample(initial_unranked_nodes, self.num_subsample))
#        elif self.rank_nodes is not None:
#            if self.verbose:
#                print("\tRunning until the given %d nodes are ranked correctly" % (len(self.rank_nodes)))
#            initial_unranked_nodes = unranked_nodes.copy()
#            # TODO if the initial set of nodes to rank is really large, we may want to use another subset
#            unranked_nodes = self.rank_nodes.copy()

        # the infinity norm is simply the maximum value in the vector
        max_f = self.f.max()
        if self.verbose:
            print("\tmax_f: %0.4f" % (max_f))

        self.num_iters = 0
        # total number of computations performed during the update function
        self.total_comp = 0
        # amount of time taken during the update function
        self.total_update_time = 0
        # TODO use cpu time, not just system time. time.process_time() should work
        start_time = time.process_time()
        # also keep track of the max score change after each iteration
        self.max_d_list = []
        # keep track of the UB after each iteration
        self.UB_list = []
        # keep track of how fast the nodes are ranked
        self.num_unranked_list = []
        self.kendalltau_list = []
        # keep track of fmax, avgp, auprc, auroc at each iteration
        self.eval_stats_list = []
        #self.spearmanr_list = []
        # keep track of the biggest # of nodes with continuously overlapping upper or lower bounds
        max_unranked_stretch = 0 
        self.max_unranked_stretch_list = [] 
        # keep track of the maximum difference of the current scores to the final score
        self.max_d_compare_ranks_list = []
        ## also keep track of how many nodes have a fixed ranking from the top of the list
        #num_ranked_from_top = 0 
        #self.num_ranked_from_top_list = [] 

        # iterate until the top-k are attained
        # R is not updated if either rank_all or rank_pos_neg is True
        while len(R) > self.k or (self.rank_topk is True and len(unranked_nodes) > 0):
            # also stop for any of these criteria
            if (self.rank_all is True or self.rank_pos_neg is not None) \
                    and len(unranked_nodes) == 0:
                # if the subset of nodes are ranked, but all nodes are not, then keep going
                if initial_unranked_nodes is not None:
                    unranked_nodes = initial_unranked_nodes
                    initial_unranked_nodes = None 
                else:
                    break
            # this was repalced with the unr_pos_nodes and unr_neg_nodes
#            elif self.rank_nodes is not None and len(unranked_nodes) == 0:
#                # if the subset is ranked, but the rank_nodes set hasn't been compared to all nodes, then keep going
#                if initial_unranked_nodes is not None:
#                    if self.verbose:
#                        print("\tGetting which nodes whose UB/LB conflict with the %d rank_nodes" % (len(self.rank_nodes)))
#                    # get the set of nodes that conflict with the rank_nodes 
#                    conflicting_nodes = alg_utils.check_fixed_rankings(LBs, UBs, initial_unranked_nodes,
#                            nodes_to_rank=self.rank_nodes) 
#                    # if there are many conflicting nodes, then just check a subset
#                    if len(conflicting_nodes) > self.num_subsample:
#                        unranked_nodes = set(random.sample(conflicting_nodes, self.num_subsample))
#                    else:
#                        conflicting_nodes = None 
#                    initial_unranked_nodes = None 
#                elif conflicting_nodes is not None:
#                    unranked_nodes = conflicting_nodes
#                    conflicting_nodes = None
#                else:
#                    break

            if self.verbose:
                print("\tnum_iters: %d, |R|: %d, |unranked_nodes|: %d, max_unranked_stretch: %d" % (
                    self.num_iters, len(R), len(unranked_nodes), max_unranked_stretch))
            if self.num_iters > self.max_iters:
                if self.verbose:
                    print("\thit the max # iters: %d. Stopping." % (self.max_iters))
                break
            self.num_iters += 1
            # keep track of how long it takes to update the bounds at each iteration
            curr_time = time.process_time()

            # power iteration
            LBs = self.a*csr_matrix.dot(self.P, prev_LBs) + self.f

            update_time = time.process_time() - curr_time
            self.total_update_time += update_time
            max_d = (LBs - prev_LBs).max()
            prev_LBs = LBs.copy()
            UB = self.computeUBs(max_f, self.a, self.num_iters)
            if self.scores_to_compare is not None:
                max_d_compare_ranks = (self.scores_to_compare - LBs).max()

            #min_delta = min([float(LBs[n] - prev_LBs[n]) for n in N])
            #print("\t\tdelta_N: %0.4f, LBs[n]: %0.4f, min_delta: %0.5f" % (delta_N, LBs[largest_diff_node], min_delta))
            if self.verbose:
                if self.scores_to_compare is not None:
                    print("\t\t%0.4f sec to update scores. max_d: %0.2e, UB: %0.2e, max_d_compare_ranks: %0.2e" % (update_time, max_d, UB, max_d_compare_ranks))
                else:
                    print("\t\t%0.4f sec to update scores. max_d: %0.2e, UB: %0.2e" % (update_time, max_d, UB))

            # check to see if the set of nodes to rank have a fixed ranking
            if self.rank_pos_neg is not None:
                UBs = LBs + UB
                self.unr_pos_nodes, self.unr_neg_nodes = alg_utils.check_fixed_rankings(
                        LBs, UBs, unr_pos_nodes=self.unr_pos_nodes, unr_neg_nodes=self.unr_neg_nodes) 
                # the sets are disjoint, so just combine them
                unranked_nodes = list(self.unr_pos_nodes) + list(self.unr_neg_nodes)

            # if the top-k nodes don't need to be ranked correctly, stop now and return the top-k
            # otherwise, continue iterating until the top-k nodes are ranked correctly
            elif (self.rank_topk is True and len(R) <= self.k) \
                 or self.rank_all is True:
                UBs = LBs + UB
                # no point checking which nodes are fixed if the UB is too high 
                if UB > 0.1:
                    pass
                    #fixed_nodes = set()
                else:
                    #fixed_nodes, max_unranked_stretch, num_ranked_from_top = alg_utils.check_fixed_rankings(LBs, UBs, unranked_nodes) 
                    unranked_nodes, max_unranked_stretch = alg_utils.check_fixed_rankings(
                            LBs, UBs, unranked_nodes=unranked_nodes) 
                #unranked_nodes = unranked_nodes - fixed_nodes

            # now check to see if there are nodes that no longer are eligible for the top-k
            elif len(R) > self.k:
                # get the score of the node with the kth largest score
                k_score = LBs[np.argpartition(LBs, -self.k)[-self.k]]

                if self.verbose:
                    print("\t\tk_score: %0.6f, additional_score: %0.6f" % (k_score, UB))

                R = [R[i] for i in np.where(LBs[R] + UB - self.epsUB >= k_score)[0]]

                if len(R) == self.k:
                    unranked_nodes = set(R.copy())

            self.max_unranked_stretch_list.append(max_unranked_stretch)
            self.max_d_list.append(max_d) 
            if self.scores_to_compare is not None:
                self.max_d_compare_ranks_list.append(max_d_compare_ranks) 
            self.UB_list.append(UB) 
            self.num_unranked_list.append(len(unranked_nodes))
            if self.ranks_to_compare is not None:
                # also append a measure of the similarity between the current ranking and the rank to compare with
                # get the current node ranks
                scores = {self.idx2node[n]:LBs[n] for n in range(len(LBs))}
                nodes_with_ranks = set(self.ranks_to_compare)
                nodes_to_rank = set(scores.keys()) & nodes_with_ranks
                # check to make sure we have a rank for all of the nodes
                if len(nodes_to_rank) != len(nodes_with_ranks):
                    print("ERROR: some nodes do not have a ranking")
                    print("\t%d nodes_to_rank, %d ranks_to_compare" % (len(nodes_to_rank), len(nodes_with_ranks)))
                    sys.exit()
                # builds a dictionary of the node as the key and the current rank as the value
                # e.g., {50: 0, 20: 1, ...}
                curr_ranks = {n: i for i, n in enumerate(sorted(nodes_to_rank, key=scores.get, reverse=True))}
                # if I sort using ranks_to_compare directly, then for the first couple iterations when many nodes are tied at 0, 
                # will be left in the order they were in (i.e., matching the correct/fixed ordering)
                #curr_ranks = {n: i for i, n in enumerate(sorted(self.ranks_to_compare, key=scores.get, reverse=True))}
                # get the current rank of the nodes in the order of the ranks_to_compare 
                # for example, if the ranks_to_compare has 20 at 0 and 50 at 1, and the current rank is 50: 0, 20: 1,
                # then compare_ranks will be [1, 0]
                compare_ranks = [curr_ranks[n] for n in self.ranks_to_compare]
                # compare the two rankings
                # for example: curr rank: [1,0], orig rank: [0,1] 
                self.kendalltau_list.append(kendalltau(compare_ranks, range(len(self.ranks_to_compare)))[0])
                # this is no longer needed
                #self.spearmanr_list.append(spearmanr(compare_ranks, range(len(self.ranks_to_compare)))[0])
                if self.rank_pos_neg is not None:
                    # need to include all nodes because otherwise the recall will be higher
                    # from the unreachable positives that were removed
                    scores_arr = np.zeros(self.num_nodes)
                    indices = [self.idx2node[n] for n in range(len(LBs))]
                    scores_arr[indices] = LBs
                    prec, recall, fpr = alg_utils.compute_eval_measures(scores_arr, self.rank_pos_neg[0], self.rank_pos_neg[1])
                    fmax = alg_utils.compute_fmax(prec, recall)
                    avgp = alg_utils.compute_avgp(prec, recall)
                    auprc = alg_utils.compute_auprc(prec, recall)
                    auroc = alg_utils.compute_auroc([r for r, f in fpr], [f for r, f in fpr])
                    self.eval_stats_list.append((fmax, avgp, auprc, auroc))

        self.total_time = time.process_time() - start_time
        self.total_comp += len(self.P.data)*self.num_iters
        #return R, LBs, total_time, num_iters, total_comp, max_d_list
        return R, LBs

    def computeUBs(self, max_f, a, i):
        if a == 1:
            return 1
        else:
            additional_score = (a**(i) * max_f) / (1-a)
        #for u in R:
        #    UBs[u] = LBs[u] + additional_score

        return additional_score

    def initialize_sets(self):
        # R is the set of candidate top-k nodes
        R = np.arange(self.P.shape[0]).astype(int)
        unranked_nodes = set(R.copy())
        #N = R.copy()

        #self.f = self.f.astype('float128')
        # set the initial lower bound (LB) of each node to f or 0
        # TODO no need to keep all nodes in the dictionary. Just the ones in B or N
        LBs = self.f.copy()
        # dictionary of LBs at the previous iteration
        #prev_LBs = self.f.copy()
        # start them at 0 so the delta_N will be the correct amount of change at each iteration
        prev_LBs = np.zeros(len(LBs))
        # dictionary of Upper Bonds for each node
        UBs = np.ones(len(R))
        #UBs = np.ones(len(R), dtype='float128')

        return R, unranked_nodes, LBs, prev_LBs, UBs

    def get_stats(self):
        """
        Returns the total time, time to update scores, # of iterations, # of computations (estimated), 
        the max_d at each iteration, and the initial size of the graph.
        """
        return self.total_time, self.total_update_time, self.num_iters, self.total_comp, self.P.shape[0]
