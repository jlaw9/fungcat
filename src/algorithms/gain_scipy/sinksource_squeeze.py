
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
from scipy.stats import kendalltau, spearmanr, weightedtau


class SinkSourceSqueeze:

    def __init__(self, P, positives, negatives=None, k=100, a=0.8, 
                 epsUB=0, rank_topk=False, rank_all=False, rank_nodes=None,
                 verbose=False, ranks_to_compare=None):
        """
        *P*: Row-normalized sparse-matrix representation of the graph
        *f*: initial vector f of amount of score received from positive nodes
        *epsUB*: if all other nodes have an UB - epsUB < the kth node's LB, then return the top-k 
        *rank_topk*: require that the ranks of the top-k nodes be fixed using their UB and LB
        *rank_all*: require that the ranks of all nodes be fixed using their UB and LB
        *rank_nodes*: set of nodes for which a fixed ranking will be required using their UB and LB.
            The k parameter will be ignored.
        *ranks_to_compare*: A dictionary of scores/ranks to compare the ranking after each iteration.
            Used to see how long it takes for the ranks to be correct, even if the ranks aren't fixed
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
        self.rank_nodes = rank_nodes
        self.ranks_to_compare = ranks_to_compare
        self.verbose = verbose

    def runSinkSourceSqueeze(self):
        # TODO this should be done once before all predictions are being made
        # check to make sure the graph is normalized because making a copy can take a long time
        #G = alg_utils.normalizeGraphEdgeWeights(G)
        self.P, self.f, self.node2idx, self.idx2node = alg_utils.setupScores(
            self.P, self.positives, self.negatives, a=self.a, 
            remove_nonreachable=True, verbose=self.verbose)
        # if rank_nodes is specified, map those node ids to the current indices
        if self.rank_nodes is not None:
            self.rank_nodes = set(self.node2idx[n] for n in self.rank_nodes if n in self.node2idx)

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
        scores = defaultdict(int)
        for n in range(len(all_LBs)):
            scores[self.idx2node[n]] = all_LBs[n]
        #scores = {self.idx2node[n]:all_LBs[n] for n in range(len(all_LBs))}

        if self.verbose or True:
            print("SinkSourceSqueeze found top k after %d iterations (%0.3f total sec, %0.3f sec to update)"
                  % (self.num_iters, self.total_time, self.total_update_time))

        return R, scores

    def _SinkSourceSqueeze(self):
        """
        *returns*: The set of top-k nodes, and current scores for all nodes
        """
        # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
        R, unranked_nodes, LBs, prev_LBs, UBs = self.initialize_sets()
        #initial_unranked_nodes = None
        if self.rank_all is True:
            if self.verbose:
                print("\tRunning until all nodes are ranked correctly")
            # see if sampling a subset of nodes will speed-up the ranking
            initial_unranked_nodes = unranked_nodes.copy()
            # randomly sample 1/100 of the nodes to check for fixed ranks
            #unranked_nodes = random.sample(initial_unranked_nodes, len(unranked_nodes) / float(100))
            unranked_nodes = set(random.sample(initial_unranked_nodes, 500))
        elif self.rank_nodes is not None:
            if self.verbose:
                print("\tRunning until the given %d nodes are ranked correctly" % (len(self.rank_nodes)))
            initial_unranked_nodes = unranked_nodes.copy()
            # TODO if the initial set of nodes to rank is really large, we may want to use another subset
            unranked_nodes = self.rank_nodes.copy()

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
        # keep track of how fast the nodes are ranked
        self.num_unranked_list = []
        self.kendalltau_list = []
        self.spearmanr_list = []

        # iterate until the top-k are attained
        while len(R) > self.k or (self.rank_topk is True and len(unranked_nodes) > 0):
            # also stop for any of these criteria
            if self.rank_all is True and len(unranked_nodes) == 0:
                # if the subset of nodes are ranked, but all nodes are not, then keep going
                if initial_unranked_nodes is not None:
                    unranked_nodes = initial_unranked_nodes
                    initial_unranked_nodes = None 
                else:
                    break
            elif self.rank_nodes is not None and len(unranked_nodes) == 0:
                # if the subset is ranked, but the rank_nodes set hasn't been compared to all nodes, then keep going
                if initial_unranked_nodes is not None:
                    if self.verbose:
                        print("\tGetting which nodes whose UB/LB conflict with the %d rank_nodes" % (len(self.rank_nodes)))
                    # get the set of nodes that conflict with the rank_nodes 
                    conflicting_nodes = alg_utils.check_fixed_rankings(LBs, UBs, initial_unranked_nodes,
                            nodes_to_rank=self.rank_nodes) 
                    # if there are many conflicting nodes, then just check a subset
                    if len(conflicting_nodes) > 500:
                        unranked_nodes = set(random.sample(conflicting_nodes, 500))
                    else:
                        conflicting_nodes = None 
                    initial_unranked_nodes = None 
                elif conflicting_nodes is not None:
                    unranked_nodes = conflicting_nodes
                    conflicting_nodes = None
                else:
                    break

            if self.verbose:
                print("\tnum_iters: %d, |R|: %d, |unranked_nodes|: %d" % (self.num_iters, len(R), len(unranked_nodes)))
            #if self.num_iters > self.max_iters:
            #    break
            self.num_iters += 1

            # keep track of how long it takes to update the bounds at each iteration
            curr_time = time.process_time()
            LBs = self.a*csr_matrix.dot(self.P, prev_LBs) + self.f
            update_time = time.process_time() - curr_time

            # TODO store this in a list and return it
            max_d = (LBs - prev_LBs).max()
            prev_LBs = LBs.copy()

            #min_delta = min([float(LBs[n] - prev_LBs[n]) for n in N])
            #print("\t\tdelta_N: %0.4f, LBs[n]: %0.4f, min_delta: %0.5f" % (delta_N, LBs[largest_diff_node], min_delta))
            self.total_update_time += update_time
            if self.verbose:
                print("\t\t%0.4f sec to update bounds. max_d: %0.4f" % (update_time, max_d))

            UB = self.computeUBs(max_f, self.a, self.num_iters)

            # check to see if the set of nodes to rank have a fixed ranking
            if self.rank_nodes is not None:
                UBs = LBs + UB
                fixed_nodes = alg_utils.check_fixed_rankings(LBs, UBs, unranked_nodes) 
                unranked_nodes = unranked_nodes - fixed_nodes

            # if the top-k nodes don't need to be ranked correctly, stop now and return the top-k
            # otherwise, continue iterating until the top-k nodes are ranked correctly
            elif (self.rank_topk is True and len(R) <= self.k) \
                 or self.rank_all is True:
                UBs = LBs + UB
                fixed_nodes = alg_utils.check_fixed_rankings(LBs, UBs, unranked_nodes) 
                unranked_nodes = unranked_nodes - fixed_nodes

            # now check to see if there are nodes that no longer are eligible for the top-k
            elif len(R) > self.k:
                # get the score of the node with the kth largest score
                k_score = LBs[np.argpartition(LBs, -self.k)[-self.k]]

                if self.verbose:
                    print("\t\tk_score: %0.6f, additional_score: %0.6f" % (k_score, UB))

                R = [R[i] for i in np.where(LBs[R] + UB - self.epsUB >= k_score)[0]]

                if len(R) == self.k:
                    unranked_nodes = set(R.copy())

            self.max_d_list.append(max_d) 
            self.num_unranked_list.append(len(unranked_nodes))
            if self.ranks_to_compare is not None:
                # also append a measure of the similarity between the current ranking and the rank to compare with
                # get the current node ranks

                scores = {self.idx2node[n]:LBs[n] for n in range(len(LBs))}
                curr_ranks = {n: i for i, n in enumerate(sorted(scores, key=scores.get, reverse=True))}
                # get the rank of the nodes in ranks_to_compare in the current ranking
                compare_ranks = [curr_ranks[n] for n in self.ranks_to_compare]
                compare_ranks = compare_ranks[:self.k] if self.rank_topk is True else compare_ranks
                #if len(self.ranks_to_compare) != len(curr_ranks):
                    # if the ranks_to_compare is a subset, then only get the ranks of those nodes(?)
                self.kendalltau_list.append(kendalltau(compare_ranks, range(len(self.ranks_to_compare)))[0])
                self.spearmanr_list.append(spearmanr(compare_ranks, range(len(self.ranks_to_compare)))[0])

        self.total_time = time.process_time() - start_time
        self.total_comp += len(self.P.data)*self.num_iters
        #return R, LBs, total_time, num_iters, total_comp, max_d_list
        return R, LBs

    def computeUBs(self, max_f, a, i):
        additional_score = (a**(i) * max_f) / (1-a)
        #for u in R:
        #    UBs[u] = LBs[u] + additional_score

        return additional_score

    def initialize_sets(self):
        # R is the set of candidate top-k nodes
        R = np.arange(self.P.shape[0]).astype(int)
        unranked_nodes = set(R.copy())
        #N = R.copy()

        # set the initial lower bound (LB) of each node to f or 0
        # TODO no need to keep all nodes in the dictionary. Just the ones in B or N
        LBs = self.f.copy()
        # dictionary of LBs at the previous iteration
        prev_LBs = self.f.copy()
        ## start them at 0 so the delta_N will be the correct amount of change at each iteration
        #prev_LBs = {n: 0 for n in R}
        # dictionary of Upper Bonds for each node
        UBs = np.ones(len(R))

        return R, unranked_nodes, LBs, prev_LBs, UBs

    def get_stats(self):
        """
        Returns the total time, time to update scores, # of iterations, # of computations (estimated), 
        the max_d at each iteration, and the initial size of the graph.
        """
        return self.total_time, self.total_update_time, self.num_iters, self.total_comp, self.P.shape[0]
