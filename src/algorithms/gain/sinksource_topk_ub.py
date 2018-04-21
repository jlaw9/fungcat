
# function to efficiently and accurately compute the top-k sinksource scores
# Algorithm and proofs adapted from:
# Zhang et al. Fast Inbound Top-K Query for Random Walk with Restart, KDD, 2015

#import SinkSource
import sys
import time
import alg_utils
import operator
# expects python 3 and networkx 2
import networkx as nx
#import numpy as np
import sinksource_ripple


def runSinkSourceTopK(G, positives, negatives=None, k=100, t=2, s=2, a=1, 
        deltaUBLB=None, epsUB=0):
    """
    I added the *a* (alpha) parameter so that I can still run SinkSourcePlus using the
    UB method of SinkSource. When *a*=1, SinkSource is not affected by it.
    """
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = normalizeGraphEdgeWeights(G)
    # TODO it takes up too much space to make a copy of the graph for each GO term and each algorithm
    H = G.copy()
    #print("Warning: deleting edges to positives in the original graph to speed up testing")
    H, f = sinksource_ripple.fixNodes(H, positives, a=a, f={}, UBs={}, LBs={})

    if negatives is not None:
        print("\t%d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, deltaUBLB=%s" \
                % (len(positives), len(negatives), k, t, s, str(a), str(deltaUBLB)))
        # also remove the negative nodes from the graph
        # assuming their scores are fixed at 0 which means they don't contribute any score to neighbors
        H.remove_nodes_from(negatives)
        for n in negatives:
            f.pop(n)
    else:
        print("\t%d positives, k=%d, t=%d, s=%d, a=%s, deltaUBLB=%s" \
                % (len(positives), k, t, s, str(a), str(deltaUBLB)))
    # Remove nodes that cannot be reached from a positive node
    # as their score will always be 0
    nonreachable_nodes = alg_utils.nonReachableNodes(H, [n for n in f if f[n] > 0])
    if len(nonreachable_nodes) > 0:
        print("\t%d nodes not reachable from a positive. Removing them from the graph" % (len(nonreachable_nodes)))
    H.remove_nodes_from(nonreachable_nodes)
    for n in nonreachable_nodes:
        f.pop(n)

#    # now relabel H to a new set of integers so I can use list indexing
#    H2, node2int, int2node = alg_utils.convert_labels_to_int(H)
#
#    # try changing f to a list
#    new_f = []
#    for n in sorted(H2.nodes()):
#        if int2node[n] in f:
#            new_f.append(f[int2node[n]])
#        else:
#            new_f.append(0)
#
#    f = new_f

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
#        R, LBs, time, iters, curr_k_score, len_N = SinkSourceTopK(H.subgraph(c), fc, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB, k_score=k_score)
#        # only update the k_score if the top-k in the connected component actually gave us k nodes
#        # many ccs are only a single node
#        if len(R) >= k:
#            k_score = curr_k_score
#        all_LBs.update(LBs)
#        #k_score = all_LBs[sorted(all_LBs, key=all_LBs.get, reverse=True)[k-1]]
#        overall_time += time
#
##    R = [int2node[n] for n in R]
##    LBs = {int2node[i]: LBs[i] for i in range(len(LBs))}
#    # to get the top-k nodes, sort the LBs which have been updated for each connected component
#    R = sorted(all_LBs, key=all_LBs.get, reverse=True)[:k]

    R, all_LBs, overall_time, iters, len_N = SinkSourceTopK(
            H, f, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB,
            epsUB=epsUB)

    return R, all_LBs, overall_time, iters, len_N


def SinkSourceTopK(G, f, k=100, t=2, s=2, a=1, deltaUBLB=None, k_score=0,
        epsUB=0):
    """
    *G*: Row-normalized Graph with only unknowns
    *f*: initial vector f of amount of score received from positive nodes
    *deltaUBLB*: cutoff for UB-LB diff to fix a node's score. If none, then node's scores will not be fixed
    *epsUB*: if all other nodes have a UB < the kth node's LB - epsUB, then return the top-k 
    *returns*: The set of top-k nodes, and current scores for all nodes
    """
    # TODO check to make sure t > 0, s > 0, k > 0, 0 < a < 1 and such
    R, N, F, B, LBs, prev_LBs, UBs = sinksource_ripple.initialize_sets(G, f)
    # Upper Bonds for each node. They start at a because that's the maximum amount they can receive from positives
    UBs = {n: a for n in R}
    prev_UBs = {n: a for n in R}

    num_iters = 0
    start_time = time.time()

    # iterate until the top-k are attained
    while len(R) > k or num_iters < 1:
        print("\tnum_iters: %d, |N|: %d, |B|: %d, |R|: %d" % (num_iters, len(N), len(B), len(R)))
        num_iters += 1

        # get the top s nodes in B with the highest LBs and UBs
        E = get_s_nodes_in_B(N, B, LBs, UBs, s)

        # Update nodes in N and B
        N, B = sinksource_ripple.update_N_B(G, N, B, E)
        F = F - N

        update_time = time.time()
        # update the scores of nodes in N t times
        for i in range(t):
            for u in N:
                LB_sum = 0
                UB_sum = 0
                for v, data in G[u].items():
                    w = data['weight']
                    LB_sum += w * prev_LBs[v]
                    UB_sum += w * prev_UBs[v]
                LBs[u] = a*LB_sum + f[u]
                UBs[u] = a*UB_sum + f[u]

            # TODO try using the current score directly
            for n in N:
                if LBs[n] < prev_LBs[n]:
                    print("Warning: LB decreased for %s: LBs[n]=%0.4f, prev_LBs[n]=%0.4f" % (n, LBs[n], prev_LBs[n]))
                if UBs[n] > prev_UBs[n]:
                    print("Warning: UB increased for %s: UBs[n]=%0.4f, prev_UBs[n]=%0.4f" % (n, UBs[n], prev_UBs[n]))
                prev_LBs[n] = LBs[n]
                prev_UBs[n] = UBs[n]
        print("\t\t%0.2f sec to update bounds" % (time.time() - update_time))

        # if there aren't at least k nodes in the vicinity of N,
        # then there is no need to prune nodes from R
        if len(N) <= k:
            continue

        # now check to see if there are nodes that no longer are eligible for the top-k
        # and if there are, remove them from R.
        # get the score of the node with the kth largest score
        R_LBs = {n:LBs[n] for n in R}
        #R_LBs = {n:LBs[n] for n in N | R}
        sorted_LBs = sorted(R_LBs, key=R_LBs.get, reverse=True)
        #print(sorted_LBs[:k+20])
        #k_score = LBs[sorted_LBs[k-1]]
        # get the score of the node with the kth largest score
        curr_k_score = LBs[sorted_LBs[k-1]]
        if curr_k_score > k_score:
            k_score = curr_k_score

        # update the UB of nodes in F to be the UB of the max boundary node
        if len(B) != 0:
            # TODO for some reason, the UB can dip below the k_score, eliminating the node from R
            # and then the node's score increases to be in the top-k. I think its node's in F whose
            # upper bound isn't set correctly
            #max_boundary_score = max([UBs[n] for n in N & B])
            max_boundary_score = max([UBs[n] for n in B])
        else:
            max_boundary_score = 0
        print("\t\tk_score: %0.3f, max_boundary_score: %0.3f" % (k_score, max_boundary_score))

        for u in F:
            # this UB will be overwritten by the next iteration using prev_UB
            UBs[u] = a*max_boundary_score
            prev_UBs[u] = a*max_boundary_score

        for u in N:
            if UBs[u] - LBs[u] < 0:
                print("Warning: UB[u] - LB[u] < 0: %0.3f. u: %s" % (UBs[u] - LBs[u], u))

        not_topk = set()
        for u in R:
            # I added the LB check to ensure none of the top-k
            # nodes are removed due to the - epsUB
            if UBs[u] - epsUB < k_score and LBs[u] < k_score:
                not_topk.add(u)
        R.difference_update(not_topk)
        #for u in R | N:
        #    if UBs[u] < k_score:
        #        R.discard(u)
        #    if UBs[u] > k_score:
        #        R.add(u)

        # check to see if all of the UBs for the kth node and up are tied
        # if so, then we have the top-k results and can quit
        # only need to check this if we're close to the end (len(R) < len(N))
        if len(R) < len(N) and len(R) > k:
            topk_plus_ties, R_UBs, k_UB = sinksource_ripple.checkTopKTies(LBs, UBs, R, k)
            if topk_plus_ties is True:
                print("\t%d node UBs are tied with the kth UB %0.5f. Finished" % (len(R_UBs), k_UB))
                break

        if deltaUBLB is not None:
            G, f, N, B, LBs, UBs = sinksource_ripple.checkFixNodes(G, f, N, B, LBs, UBs, a=a, deltaUBLB=deltaUBLB)
        #if len(R) - k < 10:
        #    print("; ".join(["%s: LB=%0.4f, UB=%0.4f" % (u, LBs[u], UBs[u]) for u in R]))

    total_time = time.time() - start_time
    print("SinkSourcePlusTopK found top k after %d iterations (%0.2f sec)" % (num_iters, total_time))

    return R, LBs, total_time, num_iters, len(N)


def get_s_nodes_in_B(N, B, LBs, UBs, s):
    """
    Get the top s nodes in B with the highest LBs and UBs
    I found that also exploring neighbors of nodes with high UBs
    speeds up convergence
    """
    if len(B) > s:
        if len(N) < s:
            B_LBs = {n:LBs[n] for n in B}
            # get the top s highest score nodes in B
            E = set(sorted(B_LBs, key=B_LBs.get, reverse=True)[:s])
        else:
            half_s = int(s/2)
            # get s nodes in B with the largest current scores
            B_LBs = {n:LBs[n] for n in B}
            # get the top s highest score nodes in B
            E = set(sorted(B_LBs, key=B_LBs.get, reverse=True)[:half_s])
            # also get the nodes in N with the highest UBs
            # to ensure that both UBs and LBs are converging
            B_UBs = {n:UBs[n] for n in B}
            E.update(set(sorted(B_UBs, key=B_UBs.get, reverse=True)[:half_s]))
    else:
        E = B.copy()

    return E
