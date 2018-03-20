
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import time
# expects python 3 and networkx 2
import networkx as nx
import numpy as np


def runSinkSourceRipple(G, positives, negatives=None, k=100, t=2, s=2, a=0.8, eps=None):
    print("\t%d positives, k=%d, t=%d, s=%d, a=%0.2f, eps=%s" % (len(positives), k, t, s, a, str(eps)))
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = normalizeGraphEdgeWeights(G)
    H = G.copy()
    #print("Warning: deleting edges to positives in the original graph to speed up testing")
    H, f = fixNodes(H, positives, a=a, f={}, UBs={}, LBs={})

    if negatives is not None:
        print("\t%d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, eps=%s" \
                % (len(positives), len(negatives), k, t, s, str(a), str(eps)))
        # also remove the negative nodes from the graph
        # assuming their scores are fixed at 0 which means they don't contribute any score to neighbors
        H.remove_nodes_from(negatives)
        for n in negatives:
            f.pop(n)
    else:
        print("\t%d positives, k=%d, t=%d, s=%d, a=%s, eps=%s" \
                % (len(positives), k, t, s, str(a), str(eps)))

    R, LBs, time, iters = SinkSourceRipple(H, f, k=k, t=t, s=s, a=a, eps=eps)

    return R, LBs, time, iters


def fixNodes(G, nodes_to_fix, a=1, f={}, UBs={}, LBs={}):
    """
    *f*: fixed amount of score/influence received from fixed nodes.
    Returns the graph with the fixed nodes removed,
    as well as the amount of score recieved from positives (assumes positives are fixed at a score of 1)
    """
    for u in G.nodes():
        if u not in f:
            f[u] = 0
    for v in nodes_to_fix:
        # for each neighbor u of v, add the score u would receive from v
        # G started out undirected, so the neighbors are the same as the predecessors
        # positives are fixed at 1
        fixed_score = 1
        # fix the nodes score at the average of the UB and LB
        if v in UBs and v in LBs:
            fixed_score = np.average((UBs[v], LBs[v]))
        for u in G.neighbors(v):
            f[u] += a*(G.edges[u,v]['weight'] * fixed_score)
            #G.remove_edge(u,v)
        G.remove_node(v)

    # now remove the fixed nodes from f, and we're finished
    for v in nodes_to_fix:
        f.pop(v)

    return G, f


#def SinkSourcePlusTopKRanked(G, f, k=100, t=2, s=2, a=0.8):
def SinkSourceRipple(G, f, k=100, t=2, s=2, a=0.8, eps=0.00001):
    """
    *G*: Row-normalized Graph with only unknowns
    *f*: initial vector f of amount of score received from positive nodes
    *eps*: cutoff for UB-LB diff to fix a node's score
    *returns*: The set of top-k nodes, and current scores for all nodes
    """
    # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
    R = set(G.nodes())
    # initialize the vicinity to be empty
    N = set()
    # initialize F to be all nodes not in the vicinity
    F = R - N
    # initialize the boundary nodes to be those with a non-zero value in f.
    # a node is a boundary node if it is in F with a score > 0,
    # or it is a node in N with a neighbor in F
    B = set([n for n in f if f[n] > 0])

    # set the initial lower bound (LB) of each node to f or 0
    # TODO no need to keep all nodes in the dictionary. Just the ones in B or N
    LBs = f.copy()
    # dictionary of LBs at the previous iteration
    prev_LBs = {n: 0 for n in R}
    # dictionary of Upper Bonds for each node
    UBs = {}

    num_iters = 0
    start_time = time.time()

    # iterate until the top-k are attained
    while len(R) > k:
        print("\tnum_iters: %d, |N|: %d, |B|: %d, |R|: %d" % (num_iters, len(N), len(B), len(R)))
        num_iters += 1
        if len(B) > s:
            # get s nodes in B with the largest current scores
            B_LBs = {n:LBs[n] for n in B}
            # get the top s highest score nodes in B
            E = sorted(B_LBs, key=B_LBs.get, reverse=True)[:s]
        else:
            E = B.copy()

        # add the nodes in E and their neighbors to N
        for u in E:
            N.add(u)
        prev_size_N = len(N)
        for u in E:
            # add u's neighbors to N
            for neighbor in G.neighbors(u):
                #neighbors_added.add(neighbor)
                N.add(neighbor)

        # update B by adding nodes in N that have have neighbors in F
        for u in N:
            if len(set(G.neighbors(u)) - N) > 0:
                B.add(u)
            elif u in B:
                B.remove(u)

        print("\t\t|E|: %d, num_neighbors added: %d" % (len(E), len(N) - prev_size_N))
        update_time = time.time()

        # update the scores of nodes in N t times
        for i in range(t):
            for u in N:
                LB_sum = 0
                # TODO only iterate over neighbors in N or B as the rest will be 0
                for v, data in G[u].items():
                    LB_sum += data['weight'] * prev_LBs[v]
                LBs[u] = a*LB_sum + f[u]

            if i == 0:
                # find the largest score difference after 1 iteration
                delta_N = 0
                largest_diff_node = None
                for n in N:
                    delta = LBs[n] - prev_LBs[n]
                    if delta > delta_N:
                        delta_N = delta
                        largest_diff_node = n
                #delta_N = max([LBs[n] - prev_LBs[n] for n in N])

            for n in N:
                prev_LBs[n] = LBs[n]
        min_delta = min([float(LBs[n] - prev_LBs[n]) for n in N])
        print("\t\tdelta_N: %0.4f, LBs[n]: %0.4f, min_delta: %0.5f" % (delta_N, LBs[largest_diff_node], min_delta))
        print("\t\t%0.2f sec to update bounds" % (time.time() - update_time))

        # if there aren't at least k nodes in the vicinity of N,
        # then there is no need to prune nodes from R
        if len(N) <= k:
            continue

        # get the score of the node with the kth largest score
        k_score = LBs[sorted(LBs, key=LBs.get, reverse=True)[k-1]]

        # now check to see if there are nodes that no longer are eligible for the top-k
        # and if there are, remove them from R
        F = F - N

        UBs = computeUBs(LBs, UBs, G, R, N, B, F, delta_N, a, t, k_score)

        not_topk = set()
        for u in R:
            if UBs[u] < k_score:
                not_topk.add(u)
        R.difference_update(not_topk)

        # check to see if all of the UBs for the kth node and up are tied
        # if so, then we have the top-k results and can quit
        # only need to check this if we're close to the end (len(R) < len(N))
        if len(R) < len(N) and len(R) > k:
            topk_plus_ties, R_UBs, k_UB = checkTopKTies(LBs, UBs, R, k)
            if topk_plus_ties is True:
                print("\t%d node UBs are tied with the kth UB %0.5f. Finished" % (len(R_UBs), k_UB))
                break

        if eps is not None:
            #checkFixNodes()
            pass
        #if len(R) - k < 10:
        #    print("; ".join(["%s: LB=%0.4f, UB=%0.4f" % (u, LBs[u], UBs[u]) for u in R]))

    total_time = time.time() - start_time
    print("SinkSourcePlusTopK found top k after %d iterations (%0.2f sec)" % (num_iters, total_time))

    return R, LBs, total_time, num_iters


def computeHopUB(max_boundary_score, a, t, delta_N, h=0):
    """
    Compute the maximum difference of score for a given node at hop *h* away from F
    TODO See X for a proof
    """
    hopUB = (a**(h+1)/(1-a**2))*max_boundary_score + \
            ((a**(h+t+1) + a**t - a**(t+2)) / (1-a-a**2+a**3))*delta_N

    if hopUB < 0:
        print("Debug: hopUB < 0: %0.2f. Setting to 0" % (hopUB))
        hopUB = 0

    #if hopUB > 1:
    #    print("Debug: hopUB > 1: %0.2f." % (hopUB))

    return hopUB


def computeUBs(LBs, UBs, G, R, N, B, F, delta_N, a, t, k_score):

    # TODO update function to only compute UBs for nodes in R
    # May not offer much of a speed-up
    # Also, if we compute the UB for each node, we can use it to see if a node's score can be fixed
    if len(B) != 0:
        max_boundary_score = max([LBs[n] for n in B])
    else:
        max_boundary_score = 0
    print("\t\tk_score: %0.3f, max_boundary_score: %0.3f" % (k_score, max_boundary_score))
    #if len(R & F) > 0:
    # first check to see if nodes in F can be removed
    F_UB = computeHopUB(max_boundary_score, a, t, delta_N, h=0)
    #if F_UB < k_score:
    #    R.difference_update(F - B)
    for u in F - B:
        UBs[u] = F_UB
    # check the nodes in F that have a LB > 0
    for u in B - N:
        UBs[u] = LBs[u] + F_UB
    # if there are nodes in F still in R,
    # then no need to check the nodes in N because they all have a LB > LBs in F
    #if len(R & F) == 0:
    #else:
    # compute the upper bound of the nodes in N using the hops
    # which is the # of steps needed to reach a node in F
    h = 1
    # start at boundary nodes in N
    curr_nodes = N & B
    curr_neighbors = set()
    nodes_left = N.copy()
    currHopUB = 0
    while len(curr_nodes) > 0:
        currHopUB = computeHopUB(max_boundary_score, a, t, delta_N, h=h)
        # continue iterating until there are no more hops
        for u in curr_nodes:
            UBs[u] = LBs[u] + currHopUB
            curr_neighbors.update(set(G.neighbors(u)) & nodes_left)
        nodes_left.difference_update(curr_nodes)
        curr_nodes = curr_neighbors
        curr_neighbors = set()
        h += 1
    if len(nodes_left) > 0:
        # TODO there could be nodes in a different connected component which are not connected to any boundary nodes
        # set their hop to infinity
        infHopUB = computeHopUB(max_boundary_score, a, t, delta_N, h=float("inf"))
        print("\t\th: %d, delta_N: %0.4f, currHopUB: %0.4f, infHopUB: %0.6f" % (h, delta_N, currHopUB, infHopUB))
        for u in nodes_left:
            UBs[u] = LBs[u] + infHopUB

    return UBs


def checkFixNodes(G, f, N, B, LBs, UBs, a=1, eps=0.00001):
    # check to see if any of the scores of nodes in N can be fixed
    # if the difference between the UB and LB is small enough
    # TODO how small does the difference need to be?
    nodes_to_fix = set()
    for u in UBs.keys() & N:
        if UBs[u] - LBs[u] < 0:
            print("Warning: UB[u] - LB[u] < 0: %0.3f" % (UBs[u] - LBs[u]))
        if UBs[u] - LBs[u] < eps:
            nodes_to_fix.add(u)
    # the fixed nodes could be top-k nodes. Keep their UB and LB if they're still in R
    if len(nodes_to_fix) > 0:
        print("%d nodes have UB - LB < %0.1e. Fixing their score and removing them from the graph" % (len(nodes_to_fix), eps))
        G, f = fixNodes(G, nodes_to_fix, a, f=f, UBs=UBs, LBs=LBs)
        # also remove the newly fixed nodes from N and B
        N.difference_update(nodes_to_fix)
        B.difference_update(nodes_to_fix)
        # also set their UB and LB to be their fixed score
        for u in nodes_to_fix:
            fixed_score = np.average((UBs[u], LBs[u]))
            UBs[u] = fixed_score
            LBs[u] = fixed_score

    return G, f, N, B, LBs, UBs


def checkTopKTies(LBs, UBs, R, k, decimal_places=4):
    # TODO make the decimal_places a parameter to the script
    R_UBs = {u:round(UBs[u], decimal_places) for u in R}
    R_LBs = {u:round(LBs[u], decimal_places) for u in R}
    sorted_k_UBs = sorted(R_UBs, key=R_UBs.get, reverse=True)[k-1:]
    sorted_k_LBs = sorted(R_LBs, key=R_LBs.get, reverse=True)[k-1:]
    #k_scores = sorted([round(LBs[u], 6) for u in R], reverse=True)[k-1:]
    if len(R) < 20:
        print(sorted([round(R_UBs[n]-R_LBs[n], 4) for n in sorted_k_UBs]))
    #print([R_LBs[n] for n in sorted_k_UBs])
    #print(k_scores)
    #print(len(R_UBs), len(sorted_k_UBs))
    k_UB = R_UBs[sorted_k_UBs[0]]
    k_LB = R_LBs[sorted_k_UBs[0]]
    topk_plus_ties = True
    # I need to check both the UB and LB to ensure the nodes are tied
    for u in sorted_k_UBs:
        if R_UBs[u] != k_UB or R_LBs[u] != k_LB:
            topk_plus_ties = False
            break
    return topk_plus_ties, R_UBs, k_UB
