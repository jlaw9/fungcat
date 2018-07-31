
# Python implementation of SinkSource
import time
from tqdm import trange, tqdm
import networkx as nx
import alg_utils
import numpy as np
import gc
from scipy.sparse import csr_matrix
#import pdb


# this is supposed to match the original implementation
def SinkSource(G, f, unknowns, max_iters=1000, delta=0.0001, a=0.8, scores={}):
    """ Iterate over the graph until all of the scores of the unknowns have converged
    *max_iters*: maximum number of iterations
    *delta*: Epsilon convergence cutoff. If all nodes scores changed by < *delta* after an iteration, then stop
    *a*: alpha converted from lambda
    *scores*: rather than start from scratch, we can start from a previous run's (or different GO term's) scores
    """
    converged = False
    iters = 0
    start_time = time.time()
    if len(scores) == 0:
        s = f.copy()
        prev_s = f.copy()
    else:
        s = scores.copy()
        not_in_network = set(G.nodes()) - set(s.keys())
        print("\t\t%d/%d passed in scores are in the network" \
                % (len(s) - len(not_in_network), len(s)))
        # set the scores for the nodes that are not in the passed in scores
        for n in not_in_network:
            if n in f:
                s[n] = f[n]
            else:
                s[n] = 0
        prev_s = s.copy()
    d = {n: 0 for n in s}

    for iters in trange(max_iters):
        for u in unknowns:
            # update the value of iterate all of the neighbors of u
            neighbor_sum = 0
            for v, vdata in G[u].items():
                neighbor_sum += vdata['weight'] * prev_s[v]

            s[u] = a*neighbor_sum + f[u]
            #data['s'] = neighbor_sum / data['degree']

        
        max_d = max(s[u] - prev_s[u] for u in unknowns)
        tqdm.write("\t\tmax score change: %0.6f" % (max_d))
        if max_d < delta:
            converged = True
#        converged = True
#        # Update the previous value, as well as how much the score changed by (delta)
#        for u in unknowns:
#            #data = G.node[u]
#            d[u] = abs(prev_s[u] - s[u])
#            # check if the score changed enough to indicate the algorithm hasn't yet converged
#            if not converged or d[u] >= delta:
#                converged = False

            #data['prev_s'] = data['s']
        for u in unknowns:
            prev_s[u] = s[u]

        #converged = hasConverged(G, unknowns, delta=delta)
        if converged:
            break

    total_time = time.time() - start_time
    print("SinkSource converged after %d iterations (%0.2f sec)" % (iters, total_time))
    #print("Scores of the top k nodes:")
    #print("\t%s" % (', '.join(["%s: %s" % (u, G.node[u]['s']) for u in list(unknowns)[:15]])))

    return s, total_time, iters


# this is supposed to match the original implementation
def SinkSource_scipy(G, f, unknowns, max_iters=1000, delta=0.0001, a=0.8, scores={}):
    """ Iterate over the graph until all of the scores of the unknowns have converged
    *max_iters*: maximum number of iterations
    *delta*: Epsilon convergence cutoff. If all nodes scores changed by < *delta* after an iteration, then stop
    *a*: alpha converted from lambda
    *scores*: rather than start from scratch, we can start from a previous run's (or different GO term's) scores
    """
    print("\tCopying G to a matrix")
    A = nx.to_scipy_sparse_matrix(G)
    del G
    gc.collect()
    #A = A.todense()
    #A = np.squeeze(np.asarray(A))
    new_f = []
    for n in sorted(f):
        new_f.append(f[n])
    f = np.asarray(new_f)
    print("\tdone")
    converged = False
    iters = 0
    if len(scores) == 0:
        s = np.asarray(f)
        prev_s = s.copy()
    else:
        s = scores.copy()
        not_in_network = set(G.nodes()) - set(s.keys())
        print("\t\t%d/%d passed in scores are in the network" \
                % (len(s) - len(not_in_network), len(s)))
        # set the scores for the nodes that are not in the passed in scores
        for n in not_in_network:
            if n in f:
                s[n] = f[n]
            else:
                s[n] = 0
        s = np.asarray(s)
        prev_s = s.copy()
    #d = {n: 0 for n in s}
    #d = np.asarray([0]*len(s))

    start_time = time.time()
    #for iters in trange(max_iters):
    for iters in range(max_iters):
        # this uses way too much ram
        #s = a*np.dot(A,prev_s) + f
        s = a*csr_matrix.dot(A,prev_s) + f
        
        max_d = (s - prev_s).max()
        print("\t\tmax score change: %0.6f" % (max_d))
        #tqdm.write("\t\tmax score change: %0.6f" % (max_d))
        if max_d < delta:
            # converged!
            break
        prev_s = s.copy()

    total_time = time.time() - start_time
    print("SinkSource converged after %d iterations (%0.2f sec)" % (iters, total_time))
    #print("Scores of the top k nodes:")
    #print("\t%s" % (', '.join(["%s: %s" % (u, G.node[u]['s']) for u in list(unknowns)[:15]])))
    # convet s back to a dictionary
    scores = {n:s[n] for n in range(len(s))}

    return scores, total_time, iters


## this is supposed to match the original implementation
#def SinkSourceUB(G, unknowns, k=10, max_iters=1000, delta=0.0001):
#    """
#    Iterate over the graph until all of the scores of the unknowns have converged
#        *max_iters*: maximum number of iterations
#    *k*: number of nodes to get scores for
#    """
#    converged = False
#    iters = 0
#    start_time = time.time()
#
#    for iters in trange(max_iters):
#        for u in unknowns:
#            data = G.node[u]
#            # update the value of iterate all of the neighbors of u
#            neighbor_sum = 0
#            # upper bound (ub) sum is computed the same way
#            ub_sum = 0
#            for v in G[u]:
#                w = G.edges[u,v]['weight']
#                neighbor_sum += w * G.node[v]['prev_s']
#                ub_sum += w * G.node[v]['prev_ub']
#
#            data['s'] = (neighbor_sum + data['c']) / data['degree']
#            data['ub'] = (ub_sum + data['c']) / data['degree']
#            #G.node[u]['s'] = neighbor_sum / float(G.degree(u, weight='weight'))
#
#        converged = True
#        # Update the previous value, as well as how much the score changed by (delta)
#        for u in unknowns:
#            data = G.node[u]
#            data['d'] = data['s'] - data['prev_s']
#            data['d_ub'] = data['prev_ub'] - data['ub']
#            # check if the score changed enough to indicate the algorithm hasn't yet converged
#            if not converged or data['d'] >= delta:
#                converged = False
#
#            data['prev_s'] = data['s']
#            data['prev_ub'] = data['ub']
#
#        #converged = hasConverged(G, unknowns, delta=delta)
#        if converged:
#            break
#
#        # see if there are any nodes to prune
#        if iters > 10:
#            # only the unknowns have the 's' and 'ub' attributes
#            scores = {n:G.node[n]['s'] for n in unknowns}
#            k_score = scores[sorted(scores, key=scores.get, reverse=True)[k]]
#            ubs = {n:G.node[n]['ub'] for n in unknowns}
#
#            nodes_to_prune = set()
#            for u in sorted(ubs, key=ubs.get):
#                if ubs[u] < k_score:
#                    nodes_to_prune.add(u)
#                else:
#                    break
#
#            if len(nodes_to_prune) > 0:
#                tqdm.write("pruning %d nodes (out of %d) at iter: %d" % (len(nodes_to_prune), len(unknowns), iters))
#                #pdb.set_trace()
#                # remove them from the unknowns, and remove their current score
#                unknowns.difference_update(nodes_to_prune)
#                for u in nodes_to_prune:
#                    data = G.node[u]
#                    # fix their score
#                    # TODO figure out the best way to fix the score
#                    data['prev_s'] = data['s'] + (data['d'] / 2.0)
#                    data['prev_ub'] = data['ub'] - (data['d_ub'] / 2.0)
#                    #del data['s']
#                    #del data['ub']
#                    # Rather than fix this node's score, which will require its neighboring nodes to continue to perform the computations,
#                    # Update the neighboring node's 'c' which is the amount of score recieved from fixed neighbors
#                    # TODO this is not setup for directed graphs
#                    for v in G[u]:
#                        w = G.edges[u,v]['weight']
#                        G.node[v]['c'] += w * data['prev_s']
#                        G.node[v]['c_ub'] += w * data['prev_ub']
#                G.remove_nodes_from(nodes_to_prune)
#
#    total_time = time.time() - start_time
#    print("SinkSource converged after %d iterations (%0.2f sec)" % (iters, total_time))
#    #print("Scores of the top k nodes:")
#    #print("\t%s" % (', '.join(["%s: %0.2f" % (u, G.node[u]['s']) for u in list(unknowns)[:k]])))
#
#    return G


#def hasConverged(G, unknowns, delta=0.0001):
#    converged = True
#    for n in unknowns:
#        d = G.node[n]['d']
#        if d >= delta:
#            converged = False
#            break
#
#    return converged


def setupScores(G, positives, negatives, unknowns, a=1):
    """
    """

    #print("Initializing scores and setting up network")
    f = {}
    # initialize all of the scores of the nodes to 0
    for u in unknowns:
        #data = G.node[u]
        #data['s'] = 0
        #data['prev_s'] = 0
        #data['ub'] = 1
        #data['prev_ub'] = 1
        #data['degree'] = float(G.degree(u, weight='weight'))
#    for n in positives:
#        #G.node[n]['s'] = 1
#        G.node[n]['prev_s'] = 1
#        #G.node[n]['ub'] = 1
#        G.node[n]['prev_ub'] = 1
#    for n in negatives:
#        #G.node[n]['s'] = 0
#        G.node[n]['prev_s'] = 0
#        #G.node[n]['ub'] = 0
#        G.node[n]['prev_ub'] = 0
#
#    for u in unknowns:
        # f is the fixed values. f_ub is the fixed amount contributing to the UB
        f[u] = 0
        #data['f_ub'] = 0
        for p in positives:
            if G.has_edge(u,p):
                f[u] += a*G.edges[u,p]['weight']
                #data['f_ub'] += G.edges[u,p]['weight']
        # don't add the negatives weight. Those only contribute to the denominator
        #for n in negatives:
        #    if G.has_edge(u,n):
        #        data['c'] += G.edges[u,n]['weight']

    # now remove the positive and negative nodes from the network
    G.remove_nodes_from(positives.union(negatives))

    return G, f


#def setupSinkSourcePlus(G, unknowns, neg_node="n", neg_weight=1):
#    """ Adds an extra negative node and connects it to all unknowns
#    with a weight of *negative_weight*
#    """
#
#    G.add_weighted_edges_from([(n, neg_node, neg_weight) for n in unknowns])
#
#    return G


#def runSinkSource(G, positives, negatives=None, k=None, max_iters=1000, delta=0.0001, a=0.8):
def runSinkSource(G, positives, negatives=None, max_iters=1000, delta=0.0001, a=0.8, scores={}):
    """
    *G*: Should already be normalized
    *positives*: set of positive nodes
    *negeatives*: set of negative nodes
    *k*: Run with upper and lower bounds to prune unnecessary nodes
    """
    # TODO this should be done once before all predictions are being made
    # check to make sure the graph is normalized because making a copy can take a long time
    #G = alg_utils.normalizeGraphEdgeWeights(G)
    H = G.copy()
    sinksourceplus = False
    if negatives is None:
        sinksourceplus = True
        #negatives = set(['n'])
        negatives = set()
    unknowns = set(G.nodes()).difference(positives).difference(negatives)

#    if sinksourceplus:
#        # set the negative node's weight to be the average of the weights connected to the positive nodes
#        neg_weight = sum([G.degree(p, weight='weight')/G.degree(p) for p in positives]) / float(len(positives))
#        print("Setting the weight of the negative node to the avg of the positive nodes: %s" % (str(neg_weight)))
#        G = setupSinkSourcePlus(G, unknowns, neg_node='n', neg_weight=neg_weight)

    H, f = setupScores(H, positives, negatives, unknowns, a=a)

    # TODO this should be done
    # replace the nodes with integers to speed things up
    #H, node2int, int2node = alg_utils.convert_labels_to_int(G)
    # the positives and negatives were already removed from the graph,
    # so everything left is an unknown
    unknowns = set(H.nodes())

#    if k is not None:
#        print("\tGetting top %d predictions" % (k))
#        G = SinkSourceUB(G, unknowns, k=k, max_iters=max_iters, delta=delta, a=a)
#    else:
    s, time, iters = SinkSource(H, f, unknowns, max_iters=max_iters, delta=delta, a=a, scores=scores)
    #s, time, iters = SinkSource_scipy(H, f, unknowns, max_iters=max_iters, delta=delta, a=a, scores=scores)

    # switch the nodes back to their names
    #G = nx.relabel_nodes(G,int2node)

    return s, time, iters
