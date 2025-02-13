
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import sys
import time
#import operator
from collections import defaultdict
# expects python 3 and networkx 2
#import networkx as nx
import numpy as np
#from profilehooks import profile
from scipy.sparse import csr_matrix
import alg_utils


class SinkSourceRipple:

    def __init__(self, P, positives, negatives=None, k=100, t=2, s=2, a=0.8, 
                 deltaUBLB=None, epsUB=0, verbose=False):
        """
        *P*: Row-normalized sparse-matrix representation of the graph
        *f*: initial vector f of amount of score received from positive nodes
        *epsUB*: if all other nodes have an UB - epsUB < the kth node's LB, then return the top-k 
        *returns*: The set of top-k nodes, and current scores for all nodes
        """
        self.P = P
        self.positives = positives
        self.negatives = negatives
        self.k = k
        self.a = a
        self.t = t
        self.s = s
        self.epsUB = epsUB
        #self.rank_topk = rank_topk
        #self.rank_nodes = rank_nodes
        self.verbose = verbose

        # a_h is a dictionary containing computations used in the upper bound
        self.a_h = {}

    def runSinkSourceRipple(self):
        # TODO this should be done once before all predictions are being made
        # check to make sure the graph is normalized because making a copy can take a long time
        #G = alg_utils.normalizeGraphEdgeWeights(G)
        self.P, self.f, self.node2idx, self.idx2node = alg_utils.setupScores(
            self.P, self.positives, self.negatives, a=self.a, 
            remove_nonreachable=True, verbose=self.verbose)
#        # if rank_nodes is specified, map those node ids to the current indices
#        if self.rank_nodes is not None:
#            self.rank_nodes = set(self.node2idx[n] for n in self.rank_nodes if n in self.node2idx)

        if self.verbose:
            if self.negatives is not None:
                print("\t%d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, epsUB=%s" % (
                    len(self.positives), len(self.negatives), self.k, self.t, self.s, str(self.a), str(self.epsUB)))
            else:
                print("\t%d positives, k=%d, t=%d, s=%d, a=%s, epsUB=%s" % (
                    len(self.positives), self.k, self.t, self.s, str(self.a), str(self.epsUB)))

        R, all_LBs = self._SinkSourceRipple()
        #R, all_LBs, overall_time, iters, total_comp, len_N, max_d_list = SinkSourceRipple(
        #        P, f, k=k, t=t, s=s, a=a, epsUB=epsUB)

        R = set([self.idx2node[n] for n in R])
        # set the default score to 0 for the rank_nodes as some of them may be unreachable from positives
        scores = defaultdict(int)
        for n in range(len(all_LBs)):
            scores[self.idx2node[n]] = all_LBs[n]
        #scores = {self.idx2node[n]:all_LBs[n] for n in range(len(all_LBs))}

        if self.verbose:
            print("SinkSourceRipple found top k after %d iterations (%0.3f total sec, %0.3f sec to update)"
                  % (self.num_iters, self.total_time, self.total_update_time))

        return R, scores

    #@profile
    def _SinkSourceRipple(self):
        """
        *returns*: The set of top-k nodes, and current scores for all nodes
        """
        # neighbors is an array of arrays containing the neighbors of each node
        #node_neighbors = alg_utils.get_neighbors(P)
        # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
        R, N, F, B, LBs, prev_LBs, UBs = self.initialize_sets()
        # easier to not use the self reference every time
        k, a, t, s = self.k, self.a, self.t, self.s
#        if self.rank_nodes is not None:
#            if self.verbose:
#                print("\tRunning until the given %d nodes are ranked correctly" % (len(self.rank_nodes)))
#            initial_unranked_nodes = unranked_nodes.copy()
#            # TODO if the initial set of nodes to rank is really large, we may want to use another subset
#            unranked_nodes = self.rank_nodes.copy()

        self.num_iters = 0
        # total number of computations performed during the update function
        self.total_comp = 0
        # amount of time taken during the update function
        self.total_update_time = 0
        start_time = time.process_time()
        # also keep track of the max score change after each iteration
        self.max_d_list = []
        # keep track of how fast the nodes are ranked
        self.num_unranked_list = []

        # iterate until the top-k are attained
        while len(R) > k:
            if self.verbose:
                print("\tnum_iters: %d, |N|: %d, |B|: %d, |R|: %d" % (self.num_iters, len(N), len(B), len(R)))
            if len(N) == 0 and len(B) == 0:
                sys.exit("Error: Size of N and B are 0")
            if len(N) < len(R) and len(B) == 0:
                print("Error: N < R (%d < %d) and B == 0" % (len(N), len(R)))
                break

            self.num_iters += 1
            E = self.get_top_s_boundary_nodes(B, LBs, s)

            # Update nodes in N and B
            # Updating using the node neighbors is slightly faster 
            # for some GO terms with a small # of annotations (like 0.05 sec),
            # but is 2-4x slower for GO terms with a large # of annotations
            #N, B, F = self.update_N_B_F(node_neighbors, N, B, E, F)
            N, B, F = self.update_N_B_F_matrix(N, B, E, F)

            curr_time = time.process_time() 
            # 2-4x faster for GO terms with a large # of annotations
            #LBs, prev_LBs, delta_N, comp = self.update_scores_full(N, F, LBs, prev_LBs, a, t)
            # 2x faster for GO terms with a small # of annotations (sinksourceplus)
            LBs, prev_LBs, delta_N, comp = self.update_scores_submatrix(N, LBs, prev_LBs, a, t)
            update_time = time.process_time() - curr_time

            if self.verbose:
                print("\t\t%0.4f sec to update bounds. delta_N: %0.2e" % (update_time, delta_N))
            self.total_comp += comp
            # it's being tracked in the update_scores function now
            # to get only the matrix multiplication time
            #self.total_update_time += update_time
            self.max_d_list.append(delta_N) 

            # if there aren't at least k nodes in the vicinity of N,
            # then there is no need to prune nodes from R
            if len(N) <= k:
                continue

            # get the score of the node with the kth largest score
            #partition = np.argpartition(LBs, -k)
            #topk = set(partition[-k:])
            #not_topk = list(sorted(R - topk))
            #k_score = LBs[partition[-k]]
            k_score = LBs[np.argpartition(LBs, -k)[-k]]

            # update the UBs
            F_UB, N_UB = self.computeUBs(LBs, UBs, 
                    R, N, B, F, delta_N, a, t, k_score)
            # TODO test using the hopsUB
            #UBs = self.computeUBs(LBs, UBs, node_neighbors, 
            #        R, N, B, F, delta_N, a, t, k_score)

            if F_UB < N_UB:
                print("\t\tF_UB: %0.4f, N_UB: %0.4f" % (F_UB, N_UB))
                continue

            N_arr = np.array(list(N))
            # now check to see if there are nodes that no longer are eligible for the top-k
            # and if there are, remove them from R
            #R = R - set(np.where(UBs[not_topk] - epsUB < k_score)[0])
            # get the set of nodes in N whose UB is >= the kth node's score
            # TODO potiential bug. Sometimes a node is removed from R when it shouldn't be
            # TODO could some nodes outside of N be in the topk?
            #R = [N_arr[i] for i in np.where(LBs_N + N_UB - epsUB >= k_score)[0]]
            R = [N_arr[i] for i in np.where(LBs[N_arr] + N_UB - self.epsUB >= k_score)[0]]

        self.total_time = time.process_time() - start_time
        # keep track of the final size of the vicinity as well
        self.len_N = len(N)

        return R, LBs
        #return R, LBs, total_time, self.num_iters, total_comp, len(N), max_d_list

    def update_scores_full(self, N, F, LBs, prev_LBs, a, t):
        Fl = np.array(list(F))

        # TODO figure out how to count the # of non-zero multiplications
        #comp = len(P_N.data) * t
        comp = 0

        # update the scores of nodes in N t times
        for i in range(t):
            # keep track of only the amount of time it takes to do the matrix multiplication
            curr_time = time.process_time()
            LBs = a*csr_matrix.dot(self.P, prev_LBs) + self.f
            self.total_update_time += time.process_time() - curr_time

            LBs[Fl] = 0
            if i == 0:
                # find the largest score difference after 1 iteration
                delta_N = (LBs - prev_LBs).max()

            prev_LBs = LBs.copy()

        return LBs, prev_LBs, delta_N, comp

    def update_scores_submatrix(self, N, LBs, prev_LBs, a, t):
        #print("\t\tGetting subnetwork")
        N_arr = np.array(list(N))
        # get the new vicinity (N) subnetwork of P
        P_N = alg_utils.select_nodes(self.P, N_arr)
        P_N = self.P[N_arr,:][:,N_arr]
        prev_LBs_N = prev_LBs[N_arr]
        f_N = self.f[N_arr]
        #print("\t\tupdating")
        comp = len(P_N.data) * t

        # update the scores of nodes in N t times
        # keep track of only the amount of time it takes to do the matrix multiplication
        for i in range(t):
            curr_time = time.process_time()
            LBs_N = a*csr_matrix.dot(P_N, prev_LBs_N) + f_N
            self.total_update_time += time.process_time() - curr_time

            if i == 0:
                # find the largest score difference after 1 iteration
                delta_N = (LBs_N - prev_LBs_N).max()

            prev_LBs_N = LBs_N.copy()
            #for n in range(len(N)):
            #    prev_LBs_N[n] = LBs_N[n]
        LBs[N_arr] = LBs_N
        prev_LBs[N_arr] = prev_LBs_N

        return LBs, prev_LBs, delta_N, comp

    def computeHopUB(self, max_boundary_score, a, t, delta_N, h=0):
        """
        Compute the maximum difference of score for a given node at hop *h* away from F
        TODO See X for a proof
        """
        # store the a and h computations to speed this up
        if (a,h) not in self.a_h:
            self.a_h[(a,h)] = (a**(h+1)/(1-a**2))
        if (a,h,t) not in self.a_h:
            self.a_h[(a,h,t)] = ((a**(h+t+1) + a**t - a**(t+2)) / (1-a-a**2+a**3))
        hopUB = self.a_h[(a,h)]*max_boundary_score + \
                self.a_h[(a,h,t)]*delta_N

        if hopUB < 0:
            print("Debug: hopUB < 0: %0.2f. Setting to 0" % (hopUB))
            hopUB = 0

        #if hopUB > 1:
        #    print("Debug: hopUB > 1: %0.2f." % (hopUB))

        return hopUB

    def computeUBs(self, LBs, UBs, R, N, B, F, delta_N, a, t, k_score):
        # May not offer much of a speed-up
        # Also, if we compute the UB for each node, we can use it to see if a node's score can be fixed
        max_boundary_score = LBs[list(B)].max() if len(B) != 0 else 0
        # TODO potential bug here. The max_boundary_score increases sometimes
        # causing the UB to increase
        #max_boundary_score = max([LBs[n] for n in N & B]) if len(B) != 0 else 0

        if self.verbose:
            print("\t\tk_score: %0.3f, max_boundary_score: %0.3f" % (k_score, max_boundary_score))
        #F_UB = 1 
        #if len(R & F) > 0:
            # first check to see if nodes in F can be removed
        F_UB = self.computeHopUB(max_boundary_score, a, t, delta_N, h=0)
            #for u in F:
            #    UBs[u] = F_UB

        # compute the upper bound of the nodes in N using the hops
        # which is the # of steps needed to reach a node in F
        # TODO compare how giving all nodes in N an upper bound with h=1 changes results
        h = 1
        N_UB = self.computeHopUB(max_boundary_score, a, t, delta_N, h=h)
        if self.verbose:
            print("\t\tN_UB: %0.4f" % (N_UB))
        #for u in R & N:
        #    UBs[u] = LBs[u] + N_UB
    #    # start at boundary nodes in N
    #    curr_nodes = B.copy()
    #    # iterate to all nodes in N until they all have their UB set
    #    nodes_left = N - B
    #    currHopUB = 0
    #    while len(curr_nodes) > 0:
    #        currHopUB = computeHopUB(max_boundary_score, a, t, delta_N, h=h)
    #
    #        # continue iterating until there are no more hops
    #        curr_neighbors = set()
    #        for u in curr_nodes:
    #            UBs[u] = LBs[u] + currHopUB
    #            curr_neighbors.update(node_neighbors[u] & nodes_left)
    #        nodes_left = nodes_left - curr_nodes
    #        curr_nodes = curr_neighbors
    #        h += 1
    #    # if there are nodes in N not reachable from boundary nodes
    #    if len(nodes_left) > 0:
    #        # set their hop to infinity
    #        infHopUB = computeHopUB(max_boundary_score, a, t, delta_N, h=float("inf"))
    #        print("\t\th: %d, delta_N: %0.4f, currHopUB: %0.4f, infHopUB: %0.6f" % (h, delta_N, currHopUB, infHopUB))
    #        nodes_left = list(nodes_left)
    #        UBs[nodes_left] = LBs[nodes_left] + infHopUB

        return F_UB, N_UB

    def get_top_s_boundary_nodes(self, B, LBs, s):
        if len(B) > s:
            # using arpartition gives a tiny improvement 
            # get the top s highest score nodes in B
            Bl = list(B)
            B_LBs = LBs[Bl]
            # argpartition is supposed to be faster than sorting
            # see here: https://stackoverflow.com/a/23734295
            # TODO get the right node ids after this operation
            E = set([Bl[i] for i in np.argpartition(B_LBs, -s)[-s:]])
            #print("%d nodes different" % (s - len(E & E2)))
        else:
            E = B.copy()
        return E

    #@profile
    # TODO this function takes ~1/2 of the total running time
    # I don't think there's anything else I can do to speed it up
    def update_N_B_F_matrix(self, N, B, E, F):
        """
        Updates B and N in place to contain the additional nodes in E
        B is updated to be the nodes in N that have neighbors in F
        """
        Fl = np.array(list(F))
        El = list(E)
        # all of E's neighbors will be added, so they will no longer be boundary nodes
        B.difference_update(E)
        prev_size_N = len(N)
        new_neighbors = set()
        # to get the neighbors of E, get the rows of E, 
        # then the nonzero columns in F 
        new_neighbors = set(Fl[self.P[El][:,Fl].getnnz(axis=0).nonzero()[0]])
        N.update(new_neighbors)
        # Update F to remove the new nodes in N
        F = F - N
        # remove false positives in the function below
        B.update(new_neighbors)
        prev_size_B = len(B)

#        # this is a little faster (.59->.39 for pathogenesis)
#        # but a little slower for GO terms with a small vicinity 
#        Nl = np.array(list(N))
#        # the neighbors of the newly added neighbors may cause some boundary nodes to no longer be boundary nodes
#        potential_non_B = B & set(Nl[self.P[list(new_neighbors)][:,Nl].getnnz(axis=0).nonzero()[0]])
#        print(len(potential_non_B))
#
#        # get the rows in N (nodes in the vicinity) and the columns in F (edges outside of N)
#        # nodes with edges outside N are the boundary nodes
#        still_boundary_nodes = self.get_boundary_nodes_matrix(np.asarray(list(potential_non_B)), list(F)) 
#        B = B - (potential_non_B - still_boundary_nodes)

        # get the rows in N (nodes in the vicinity) and the columns in F (edges outside of N)
        # nodes with edges outside N are the boundary nodes
        B = self.get_boundary_nodes_matrix(np.asarray(list(B)), list(F))

        if self.verbose is True:
            print("\t\t|E|: %d, num_neighbors added: %d, num B removed: %d" % (
                len(E), len(N) - prev_size_N, prev_size_B - len(B)))

        return N, B, F

    # TODO figure out which one is faster
    def update_N_B_F(self, node_neighbors, N, B, E, F, verbose=True):
        """
        Updates B and N in place to contain the additional nodes in E
        B is updated to be the nodes in N that have neighbors in F
        """
        # all of E's neighbors will be added, so they will no longer be boundary nodes
        B.difference_update(E)
        prev_size_N = len(N)
        new_neighbors = set()
        for u in E:
            new_neighbors.update(node_neighbors[u]) 
            # add u's neighbors to N
            #N.update(node_neighbors[u])

        N.update(new_neighbors)
        # remove false positives in the loop below
        B.update(new_neighbors)
        # Update F to remove the new nodes in N
        F = F - N

        prev_size_B = len(B)
        # loop through the boundary nodes and newly added nodes to see 
        # which of them are boundary nodes
        for u in B.copy():
            boundary_node = False
            for i in node_neighbors[u]:
                if i not in N:
                    boundary_node = True 
                    B.add(u)
                    break
            if not boundary_node: 
                B.discard(u)
    #    B = get_boundary_nodes(node_neighbors, N)
        if verbose is True:
            print("\t\t|E|: %d, num_neighbors added: %d, num B removed: %d" % (
                len(E), len(N) - prev_size_N, prev_size_B - len(B)))

        return N, B, F

    def get_boundary_nodes(self, node_neighbors, N):
        #print("\t\tGetting boundary nodes")
        #B = np.array([]).astype(int)
        B = []
        # update B by adding nodes in N that have have neighbors in F
        # We need to check each node because there could be nodes in B
        # that are no longer boundary nodes from neighbors of nodes in E being added
        # TODO this should be faster, but for some reason it isn't working
        #for u in neighbors_added | B:
        for u in N:
            #neighbors = node_neighbors[u]
            #if len(np.setdiff1d(neighbors, N, assume_unique=True)) > 0:
            if len(node_neighbors[u] - N) > 0:
                B.append(u)

        return set(B)

    def get_boundary_nodes_matrix(self, N, F):
        #print("\t\tGetting boundary nodes")

        B = set(N[self.P[N][:,F].getnnz(axis=1).nonzero()[0]])

        return B

    def initialize_sets(self):
        """
        Initializes all of the node sets and score dictionaries
        """
        # R is the set of candidate top-k nodes
        # python set operations are faster than numpy 
        R = set(range(self.P.shape[0]))
        # Initialize the vicinity to be the nodes with a non-zero score
        # Otherwise if a node in B (non-zero score) that's not in N has the maximum boundary score,
        # it's score could increase after the next iteration causing the upper bound to increase
        N = set(np.nonzero(self.f)[0].astype(int))
        # initialize F to be all nodes not in the vicinity
        F = R - N
        # Boundary nodes are nodes in N with a neighbor in F
        # TODO this is too slow...
        B = self.get_boundary_nodes_matrix(np.asarray(list(N)), list(F))
        #B = np.setdiff1d(csr_matrix.dot(P, f).nonzero()[0], f.nonzero()[0])

        # set the initial lower bound (LB) of each node to f or 0
        # TODO no need to keep all nodes in the array. Just the ones in B or N
        LBs = self.f.copy()
        # dictionary of LBs at the previous iteration
        # start them at 0 so the delta_N will be the correct amount of change at each iteration
        prev_LBs = np.zeros(len(R))

        # dictionary of Upper Bonds for each node
        UBs = np.ones(len(R))

        return R, N, F, B, LBs, prev_LBs, UBs

    def get_stats(self):
        """
        Returns the total time, time to update scores, # of iterations, # of computations (estimated), 
        the max_d at each iteration, and the size of the vicinity.
        """
        return self.total_time, self.total_update_time, self.num_iters, self.total_comp, self.len_N
