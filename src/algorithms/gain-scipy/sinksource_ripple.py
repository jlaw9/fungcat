
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import sys
import time
import operator
# expects python 3 and networkx 2
#import networkx as nx
import numpy as np
#from profilehooks import profile
from scipy.sparse import csr_matrix
import alg_utils


def runSinkSourceRipple(P, positives, negatives=None, k=100, t=2, s=2, a=0.8, 
        deltaUBLB=None, epsUB=0):
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = alg_utils.normalizeGraphEdgeWeights(G)
    P, f, int2int = alg_utils.setupScores(P, positives, negatives, a=a)

    if negatives is not None:
        print("\t%d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, epsUB=%s" \
                % (len(positives), len(negatives), k, t, s, str(a), str(epsUB)))
    else:
        print("\t%d positives, k=%d, t=%d, s=%d, a=%s, epsUB=%s" \
                % (len(positives), k, t, s, str(a), str(epsUB)))
#    # Remove nodes that cannot be reached from a positive node
#    # as their score will always be 0
#    nonreachable_nodes = alg_utils.nonReachableNodes(H, [n for n in f if f[n] > 0])
#    if len(nonreachable_nodes) > 0:
#        print("\t%d nodes not reachable from a positive. Removing them from the graph" % (len(nonreachable_nodes)))
#    H.remove_nodes_from(nonreachable_nodes)
#    for n in nonreachable_nodes:
#        f.pop(n)

    R, all_LBs, overall_time, iters, total_comp, len_N, max_d_list = SinkSourceRipple(
            P, f, k=k, t=t, s=s, a=a, epsUB=epsUB)

    R = set([int2int[n] for n in R])
    all_LBs = {int2int[n]:all_LBs[n] for n in range(len(all_LBs))}

    return R, all_LBs, overall_time, iters, total_comp, len_N, max_d_list


#@profile
def SinkSourceRipple(P, f, k=100, t=2, s=2, a=0.8, epsUB=0):
    """
    *P*: Row-normalized sparse-matrix representation of the graph
    *f*: initial vector f of amount of score received from positive nodes
    *deltaUBLB*: cutoff for UB-LB diff to fix a node's score
    *epsUB*: if all other nodes have an UB - epsUB < the kth node's LB, then return the top-k 
    *returns*: The set of top-k nodes, and current scores for all nodes
    """
    # neighbors is an array of arrays containing the neighbors of each node
    #node_neighbors = alg_utils.get_neighbors(P)
    # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
    R, N, F, B, LBs, prev_LBs, UBs = initialize_sets(P, f)

    num_iters = 0
    # total number of computations performed during the update function
    total_comp = 0
    start_time = time.time()
    max_d_list = []

    # iterate until the top-k are attained
    while len(R) > k:
        print("\tnum_iters: %d, |N|: %d, |B|: %d, |R|: %d" % (num_iters, len(N), len(B), len(R)))
        if len(N) == 0 and len(B) == 0:
            sys.exit("Error: Size of N and B are 0")
        if len(N) < len(R) and len(B) == 0:
            print("Error: N < R (%d < %d) and B == 0" % (len(N), len(R)))
            break

        num_iters += 1
        E = get_top_s_boundary_nodes(B, LBs, s)

        # Update nodes in N and B
        # Updating using the node neighbors is slightly faster 
        # for some GO terms with a small # of annotations (like 0.05 sec),
        # but is 2-4x slower for GO terms with a large # of annotations
        #N, B, F = update_N_B_F(P, node_neighbors, N, B, E, F)
        N, B, F = update_N_B_F_matrix(P, N, B, E, F)

        update_time = time.time()
        # 2-4x faster for GO terms with a large # of annotations
        #LBs, prev_LBs, delta_N, comp = update_scores_full(P, f, N, F, LBs, prev_LBs, a, t)
        # 2x faster for GO terms with a small # of annotations (sinksourceplus)
        LBs, prev_LBs, delta_N, comp = update_scores_submatrix(P, f, N, LBs, prev_LBs, a, t)
        total_comp += comp

        print("\t\t%0.2f sec to update bounds. delta_N: %0.4f" % (time.time() - update_time, delta_N))
        max_d_list.append(delta_N) 

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
        F_UB, N_UB = computeUBs(LBs, UBs, 
                R, N, B, F, delta_N, a, t, k_score)
        # TODO test using the hopsUB
        #UBs = computeUBs(LBs, UBs, node_neighbors, 
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
        R = [N_arr[i] for i in np.where(LBs[N_arr] + N_UB - epsUB >= k_score)[0]]

    total_time = time.time() - start_time
    print("SinkSourcePlusTopK found top k after %d iterations (%0.2f sec)" % (num_iters, total_time))

    return R, LBs, total_time, num_iters, total_comp, len(N), max_d_list


def update_scores_full(P, f, N, F, LBs, prev_LBs, a, t):
    Fl = np.array(list(F))

    #print("\t\tGetting subnetwork")
    # get the new vicinity (N) subnetwork of P
    #P_N = alg_utils.select_nodes(P, N_arr)
    #P_N = P[N_arr,:][:,N_arr]
    #prev_LBs_N = prev_LBs[N_arr]
    #prev_LBs_N = prev_LBs.copy()
    #prev_LBs[list(F)] = 0
    ##LBs_N = LBs[N]
    #f_N = f[N_arr]
    f_N = f
    #print("\t\tupdating")
    # TODO figure out how to count the # of non-zero multiplications
    #comp = len(P_N.data) * t
    comp = 0

    # update the scores of nodes in N t times
    for i in range(t):
        #LBs_N = a*csr_matrix.dot(P_N, prev_LBs_N) + f_N
        LBs = a*csr_matrix.dot(P, prev_LBs) + f_N

        LBs[Fl] = 0
        if i == 0:
            # find the largest score difference after 1 iteration
            delta_N = (LBs - prev_LBs).max()

        prev_LBs = LBs.copy()
        #prev_LBs_N = LBs_N.copy()

    return LBs, prev_LBs, delta_N, comp


def update_scores_submatrix(P, f, N, LBs, prev_LBs, a, t):
    #print("\t\tGetting subnetwork")
    N_arr = np.array(list(N))
    # get the new vicinity (N) subnetwork of P
    P_N = alg_utils.select_nodes(P, N_arr)
    P_N = P[N_arr,:][:,N_arr]
    prev_LBs_N = prev_LBs[N_arr]
    f_N = f[N_arr]
    #print("\t\tupdating")
    comp = len(P_N.data) * t

    # update the scores of nodes in N t times
    for i in range(t):
        LBs_N = a*csr_matrix.dot(P_N, prev_LBs_N) + f_N

        if i == 0:
            # find the largest score difference after 1 iteration
            delta_N = (LBs_N - prev_LBs_N).max()

        prev_LBs_N = LBs_N.copy()
        #for n in range(len(N)):
        #    prev_LBs_N[n] = LBs_N[n]
    LBs[N_arr] = LBs_N
    prev_LBs[N_arr] = prev_LBs_N

    return LBs, prev_LBs, delta_N, comp


def computeHopUB(max_boundary_score, a, t, delta_N, h=0):
    """
    Compute the maximum difference of score for a given node at hop *h* away from F
    TODO See X for a proof
    """
    global a_h
    # store the a and h computations to speed this up
    if (a,h) not in a_h:
        a_h[(a,h)] = (a**(h+1)/(1-a**2))
    if (a,h,t) not in a_h:
        a_h[(a,h,t)] = ((a**(h+t+1) + a**t - a**(t+2)) / (1-a-a**2+a**3))
    hopUB = a_h[(a,h)]*max_boundary_score + \
            a_h[(a,h,t)]*delta_N

    if hopUB < 0:
        print("Debug: hopUB < 0: %0.2f. Setting to 0" % (hopUB))
        hopUB = 0

    #if hopUB > 1:
    #    print("Debug: hopUB > 1: %0.2f." % (hopUB))

    return hopUB


def computeUBs(LBs, UBs, R, N, B, F, delta_N, a, t, k_score):
    # May not offer much of a speed-up
    # Also, if we compute the UB for each node, we can use it to see if a node's score can be fixed
    max_boundary_score = LBs[list(B)].max() if len(B) != 0 else 0
    # TODO potential bug here. The max_boundary_score increases sometimes
    # causing the UB to increase
    #max_boundary_score = max([LBs[n] for n in N & B]) if len(B) != 0 else 0

    print("\t\tk_score: %0.3f, max_boundary_score: %0.3f" % (k_score, max_boundary_score))
    #F_UB = 1 
    #if len(R & F) > 0:
        # first check to see if nodes in F can be removed
    F_UB = computeHopUB(max_boundary_score, a, t, delta_N, h=0)
        #for u in F:
        #    UBs[u] = F_UB

    # compute the upper bound of the nodes in N using the hops
    # which is the # of steps needed to reach a node in F
    # TODO compare how giving all nodes in N an upper bound with h=1 changes results
    h = 1
    N_UB = computeHopUB(max_boundary_score, a, t, delta_N, h=h)
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


def get_top_s_boundary_nodes(B, LBs, s):
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
def update_N_B_F_matrix(P, N, B, E, F, verbose=True):
    """
    Updates B and N in place to contain the additional nodes in E
    B is updated to be the nodes in N that have neighbors in F
    """
    Fl = list(F)
    El = list(E)
    # all of E's neighbors will be added, so they will no longer be boundary nodes
    B.difference_update(E)
    prev_size_N = len(N)
    new_neighbors = set()
    # to get the neighbors of E, get the rows of E, 
    # then the nonzero columns in F 
    new_neighbors = set([Fl[i] for i in P[El][:,Fl].sum(axis=0).nonzero()[1]])
    N.update(new_neighbors)
    # Update F to remove the new nodes in N
    F = F - N
    # remove false positives in the loop below
    B.update(new_neighbors)

    prev_size_B = len(B)
    # get the rows in N (nodes in the vicinity) and the columns in F (edges outside of N)
    # nodes with edges outside N are the boundary nodes
    B = get_boundary_nodes_matrix(P, list(B), list(F))

    if verbose is True:
        print("\t\t|E|: %d, num_neighbors added: %d, num B removed: %d" % (
            len(E), len(N) - prev_size_N, prev_size_B - len(B)))

    return N, B, F


# TODO figure out which one is faster
def update_N_B_F(P, node_neighbors, N, B, E, F, verbose=True):
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

def get_boundary_nodes(node_neighbors, N):
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


def get_boundary_nodes_matrix(P, N, F):
    #print("\t\tGetting boundary nodes")

    B = set([N[i] for i in P[N][:,F].sum(axis=1).nonzero()[0]])

    return B


def initialize_sets(P, f):
    """
    Initializes all of the node sets and score dictionaries
    """
    # a_h is a dictionary containing computations used in the upper bound
    global a_h
    a_h = {}

    # R is the set of candidate top-k nodes
    # python set operations are faster than numpy 
    R = set(range(P.shape[0]))
    # Initialize the vicinity to be the nodes with a non-zero score
    # Otherwise if a node in B (non-zero score) that's not in N has the maximum boundary score,
    # it's score could increase after the next iteration causing the upper bound to increase
    N = set(np.nonzero(f)[0].astype(int))
    # initialize F to be all nodes not in the vicinity
    F = R - N
    # Boundary nodes are nodes in N with a neighbor in F
    # TODO this is too slow...
    B = get_boundary_nodes_matrix(P, list(N), list(F))
    #B = np.setdiff1d(csr_matrix.dot(P, f).nonzero()[0], f.nonzero()[0])

    # set the initial lower bound (LB) of each node to f or 0
    # TODO no need to keep all nodes in the array. Just the ones in B or N
    LBs = f.copy()
    # dictionary of LBs at the previous iteration
    # start them at 0 so the delta_N will be the correct amount of change at each iteration
    prev_LBs = np.zeros(len(R))

    # dictionary of Upper Bonds for each node
    UBs = np.ones(len(R))

    return R, N, F, B, LBs, prev_LBs, UBs

