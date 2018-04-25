
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import sys
import time
import operator
# expects python 3 and networkx 2
import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix
import alg_utils


def runSinkSourceSqueeze(P, positives, negatives=None, k=100, a=0.8, 
        epsUB=0, ranked=False, verbose=False):
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = alg_utils.normalizeGraphEdgeWeights(G)
    P, f, int2int = alg_utils.setupScores(P, positives, negatives, a=a)

    if verbose:
        if negatives is not None:
            print("\t%d positives, %d negatives, %d unknowns, k=%d, a=%s, epsUB=%s" \
                    % (len(positives), len(negatives), P.shape[0], k, str(a), str(epsUB)))
        else:
            print("\t%d positives, %d unknowns, k=%d, a=%s, epsUB=%s" \
                    % (len(positives), P.shape[0], k, str(a), str(epsUB)))

    R, all_LBs, overall_time, iters, comp, max_ds = SinkSourceSqueeze(
            P, f, k=k, a=a, epsUB=epsUB, ranked=ranked)

    R = set([int2int[n] for n in R])
    all_LBs = {int2int[n]:all_LBs[n] for n in range(len(all_LBs))}

    if verbose:
        print("SinkSourceSqueeze found top k after %d iterations (%0.2f sec)" % (iters, overall_time))

    return R, all_LBs, overall_time, iters, comp, max_ds


def SinkSourceSqueeze(P, f, max_iters=1000, k=100, a=0.8, epsUB=0, ranked=False, verbose=False):
    """
    *G*: Row-normalized Graph with only unknowns
    *f*: initial vector f of amount of score received from positive nodes
    *deltaUBLB*: cutoff for UB-LB diff to fix a node's score
    *ranked*: require that the ranks of the top-k nodes be fixed using their UB and LB
    *returns*: The set of top-k nodes, and current scores for all nodes
    """
    # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
    R, unranked_nodes, LBs, prev_LBs, UBs = initialize_sets(P, f)

    # the infinity norm is simply the maximum value in the vector
    max_f = f.max()
    if verbose:
        print("\tmax_f: %0.4f" % (max_f))

    # total number of computations performed during the update function
    total_comp = 0
    num_iters = 0
    start_time = time.time()
    # also keep track of the max score change after each iteration
    max_d_list = []

    # iterate until the top-k are attained
    while len(R) > k or (ranked is True and len(unranked_nodes) > 0):
        if verbose:
            print("\tnum_iters: %d, |R|: %d, |unranked_nodes|: %d" % (num_iters, len(R), len(unranked_nodes)))
        if num_iters > max_iters:
            break
        num_iters += 1

        update_time = time.time()
        LBs = a*csr_matrix.dot(P, prev_LBs) + f

        # TODO store this in a list and return it
        max_d = (LBs - prev_LBs).max()
        max_d_list.append(max_d) 
        prev_LBs = LBs.copy()

        #min_delta = min([float(LBs[n] - prev_LBs[n]) for n in N])
        #print("\t\tdelta_N: %0.4f, LBs[n]: %0.4f, min_delta: %0.5f" % (delta_N, LBs[largest_diff_node], min_delta))
        if verbose:
            print("\t\t%0.2f sec to update bounds. max_d: %0.4f" % (time.time() - update_time, max_d))

        # get the score of the node with the kth largest score
        k_score = LBs[np.argpartition(LBs, -k)[-k]]

        # now check to see if there are nodes that no longer are eligible for the top-k
        if len(R) > k:
            UB = computeUBs(LBs, UBs, P, R, max_f, a, num_iters, k_score)
            if verbose:
                print("\t\tk_score: %0.6f, additional_score: %0.6f" % (k_score, UB))

            R = [R[i] for i in np.where(LBs[R] + UB - epsUB >= k_score)[0]]

            if len(R) == k:
                unranked_nodes = R.copy()

        # if the top-k nodes don't need to be ranked correctly, stop now and return the top-k
        # otherwise, continue iterating until the top-k nodes are ranked correctly
        if ranked is True and len(unranked_nodes) <= k:
            UBs = computeUBs(LBs, UBs, G, unranked_nodes, max_f, a, num_iters, k_score)
            # find all of the nodes in the top-k whose rankings are fixed
            fixed_ranks = {u: True for u in unranked_nodes}
            for u in unranked_nodes:
                for v in unranked_nodes:
                    # don't check the same pair of nodes twice
                    if v <= u:
                        continue
                    # if the LB and UB of these two nodes overlap, then these two node's rankings are not fixed
                    if LBs[u] < UBs[v] and UBs[u] > LBs[v]:
                        fixed_ranks[u] = False
                        fixed_ranks[v] = False
                        if len(unranked_nodes) < 70:
                            print("\t\tLBs[u] < UBs[v] and UBs[u] > LBs[v]: %0.6f < %0.6f, %0.6f > %0.6f" % (LBs[u], UBs[v], UBs[u], LBs[v]))
                        break
            for u in fixed_ranks:
                if fixed_ranks[u] is True:
                    unranked_nodes.discard(u)

            if len(unranked_nodes) < 70:
                print("\t\tnodes whose ranks aren't fixed:")
                sorted_u = sorted(unranked_nodes, key=LBs.get, reverse=True)
                print("; ".join(["UB:%0.6f,LB:%0.6f" % (UBs[n], LBs[n]) for n in unranked_nodes]))
            else:
                print("\t\t%d node ranks aren't fixed" % (len(unranked_nodes)))

    total_time = time.time() - start_time
    total_comp += len(P.data)*num_iters
    return R, LBs, total_time, num_iters, total_comp, max_d_list


def computeUBs(LBs, UBs, P, R, max_f, a, i, k_score):
    additional_score = (a**(i) * max_f) / (1-a)
    #for u in R:
    #    UBs[u] = LBs[u] + additional_score

    return additional_score


def initialize_sets(P, f):
    # R is the set of candidate top-k nodes
    R = np.arange(P.shape[0]).astype(int)
    unranked_nodes = R.copy()
    #N = R.copy()

    # set the initial lower bound (LB) of each node to f or 0
    # TODO no need to keep all nodes in the dictionary. Just the ones in B or N
    LBs = f.copy()
    # dictionary of LBs at the previous iteration
    prev_LBs = f.copy()
    ## start them at 0 so the delta_N will be the correct amount of change at each iteration
    #prev_LBs = {n: 0 for n in R}
    # dictionary of Upper Bonds for each node
    UBs = np.ones(len(R))

    return R, unranked_nodes, LBs, prev_LBs, UBs
