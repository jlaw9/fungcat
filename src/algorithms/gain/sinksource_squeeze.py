
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
import alg_utils


def runSinkSourceSqueeze(G, positives, negatives=None, k=100, a=0.8, deltaUBLB=None, ranked=False):
    #print("\t%d positives, k=%d, t=%d, s=%d, a=%0.2f, deltaUBLB=%s" % (len(positives), k, t, s, a, str(deltaUBLB)))
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = alg_utils.normalizeGraphEdgeWeights(G)
    H = G.copy()
    #print("Warning: deleting edges to positives in the original graph to speed up testing")
    H, f = fixNodes(H, positives, a=a, f={}, UBs={}, LBs={})

    if negatives is not None:
        # also remove the negative nodes from the graph
        # assuming their scores are fixed at 0 which means they don't contribute any score to neighbors
        H.remove_nodes_from(negatives)
        for n in negatives:
            f.pop(n)
        #print("\t%d positives, %d negatives, %d unknowns, k=%d, a=%s, deltaUBLB=%s" \
        #        % (len(positives), len(negatives), H.number_of_nodes(), k, str(a), str(deltaUBLB)))
        print("\t%d positives, %d negatives, %d unknowns, k=%d, a=%s" \
                % (len(positives), len(negatives), H.number_of_nodes(), k, str(a)))
    else:
        #print("\t%d positives, %d unknowns, k=%d, a=%s, deltaUBLB=%s" \
        #        % (len(positives), H.number_of_nodes(), k, str(a), str(deltaUBLB)))
        print("\t%d positives, %d unknowns, k=%d, a=%s" \
                % (len(positives), H.number_of_nodes(), k, str(a)))
    # Remove nodes that cannot be reached from a positive node
    # as their score will always be 0
    nonreachable_nodes = alg_utils.nonReachableNodes(H, [n for n in f if f[n] > 0])
    if len(nonreachable_nodes) > 0:
        print("\t%d nodes not reachable from a positive. Removing them from the graph" % (len(nonreachable_nodes)))
    H.remove_nodes_from(nonreachable_nodes)
    for n in nonreachable_nodes:
        f.pop(n)

    # Section for computing topk on each connected component
#    overall_time = 0
#    all_LBs = {}
#    k_score = 0
#    # now split the graph into connected components (ccs), and find the top-k for each component.
#    # TODO If the UB of the nodes in the new connected component are < LB for the top-k from previous CC, then stop
#    ccs = list(nx.weakly_connected_components(H))
#    print("%d connected components in H" % (len(ccs)))
#    # TODO 
#    # find the cc with the node with the largest influence from positive nodes, and start there
#    while len(ccs) > 0:
#        max_f = max(f.items(), key=operator.itemgetter(1))[0]
#        for i, c in enumerate(ccs):
#            if max_f in c:
#                break
#    #for c in ccs:
#        print("Size of connected component with largest f value: %d " % (len(c)))
#        
#        del ccs[i]
#        fc = {n:f[n] for n in f.keys() & c}
#        for n in fc:
#            f.pop(n)
#        R, LBs, time, iters, k_score = SinkSourceSqueeze(H.subgraph(c), fc, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB, k_score=k_score)
#        all_LBs.update(LBs) 
#        #k_score = all_LBs[sorted(all_LBs, key=all_LBs.get, reverse=True)[k-1]]
#        overall_time += time
#    #R = sorted(all_LBs, key=all_LBs.get, reverse=True)[:k]
#    overall_time += time

    R, all_LBs, overall_time, iters, len_N = SinkSourceSqueeze(H, f, k=k, a=a, ranked=ranked)

    return R, all_LBs, overall_time, iters, len_N


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
def SinkSourceSqueeze(G, f, k=100, a=0.8, ranked=False):
    """
    *G*: Row-normalized Graph with only unknowns
    *f*: initial vector f of amount of score received from positive nodes
    *deltaUBLB*: cutoff for UB-LB diff to fix a node's score
    *ranked*: require that the ranks of the top-k nodes be fixed using their UB and LB
    *returns*: The set of top-k nodes, and current scores for all nodes
    """
    # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
    R = set(G.nodes())
    unranked_nodes = set(G.nodes())
    N = set(G.nodes())

    # the infinity norm is simply the maximum value in the vector
    max_f = max(f.items(), key=operator.itemgetter(1))[1]
    print("\tmax_f: %0.4f" % (max_f))

    # set the initial lower bound (LB) of each node to f or 0
    # TODO no need to keep all nodes in the dictionary. Just the ones in B or N
    LBs = f.copy()
    # dictionary of LBs at the previous iteration
    prev_LBs = f.copy()
    ## start them at 0 so the delta_N will be the correct amount of change at each iteration
    #prev_LBs = {n: 0 for n in R}
    # dictionary of Upper Bonds for each node
    UBs = {}

    num_iters = 0
    start_time = time.time()

    # iterate until the top-k are attained
    while len(R) > k or (ranked is True and len(unranked_nodes) > 0):
        print("\tnum_iters: %d, |R|: %d, |unranked_nodes|: %d" % (num_iters, len(R), len(unranked_nodes)))
        #if len(N) == 0 and len(B) == 0:
        #    sys.exit("Error: Size of N and B are 0")

        num_iters += 1
        update_time = time.time()

        for u in N:
            LB_sum = 0
            # TODO only iterate over neighbors in N or B as the rest will be 0
            for v, data in G[u].items():
                LB_sum += data['weight'] * prev_LBs[v]
            LBs[u] = a*LB_sum + f[u]

        for n in N:
            prev_LBs[n] = LBs[n]
        #min_delta = min([float(LBs[n] - prev_LBs[n]) for n in N])
        #print("\t\tdelta_N: %0.4f, LBs[n]: %0.4f, min_delta: %0.5f" % (delta_N, LBs[largest_diff_node], min_delta))
        print("\t\t%0.2f sec to update bounds" % (time.time() - update_time))

        # get the score of the node with the kth largest score
        sorted_LBs = sorted(LBs, key=LBs.get, reverse=True)
        k_score = LBs[sorted_LBs[k-1]]

        # now check to see if there are nodes that no longer are eligible for the top-k
        if len(R) > k:
            UBs = computeUBs(LBs, UBs, G, R, max_f, a, num_iters, k_score)
            not_topk = set()
            for u in R:
                if UBs[u] < k_score:
                    not_topk.add(u)
            R.difference_update(not_topk)

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

#        # check to see if all of the UBs for the kth node and up are tied
#        # if so, then we have the top-k results and can quit
#        # only need to check this if we're close to the end (len(R) < len(N))
#        if len(R) < len(N) and len(R) > k:
#            topk_plus_ties, R_UBs, k_UB = checkTopKTies(LBs, UBs, R, k)
#            if topk_plus_ties is True:
#                print("\t%d node UBs are tied with the kth UB %0.5f. Finished" % (len(R_UBs), k_UB))
#                break

        #if deltaUBLB is not None:
        #    G, f, N, B, LBs, UBs = checkFixNodes(G, f, N, B, LBs, UBs, a=a, deltaUBLB=deltaUBLB)
        #if len(R) - k < 10:
        #    print("; ".join(["%s: LB=%0.4f, UB=%0.4f" % (u, LBs[u], UBs[u]) for u in R]))

    total_time = time.time() - start_time
    print("SinkSourcePlusTopK found top k after %d iterations (%0.2f sec)" % (num_iters, total_time))

    return R, LBs, total_time, num_iters, len(N)


def computeUBs(LBs, UBs, G, R, max_f, a, i, k_score):
    additional_score = (a**(i) * max_f) / (1-a)
    print("\t\tk_score: %0.6f, additional_score: %0.6f" % (k_score, additional_score))
    for u in R:
        UBs[u] = LBs[u] + additional_score

    return UBs


def checkFixNodes(G, f, N, B, LBs, UBs, a=1, deltaUBLB=0.00001):
    # check to see if any of the scores of nodes in N can be fixed
    # if the difference between the UB and LB is small enough
    # TODO how small does the difference need to be?
    nodes_to_fix = set()
    for u in UBs.keys() & N:
        if UBs[u] - LBs[u] < 0:
            print("Warning: UB[u] - LB[u] < 0: %0.3f" % (UBs[u] - LBs[u]))
        if UBs[u] - LBs[u] < deltaUBLB:
            nodes_to_fix.add(u)
    # the fixed nodes could be top-k nodes. Keep their UB and LB if they're still in R
    if len(nodes_to_fix) > 0:
        print("%d nodes have UB - LB < %0.1e. Fixing their score and removing them from the graph" % (len(nodes_to_fix), deltaUBLB))
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
