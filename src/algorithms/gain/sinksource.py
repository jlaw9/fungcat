
# Python implementation of SinkSource

#import networkx as nx
import time
from tqdm import trange, tqdm
import networkx as nx
#import pdb


# this is supposed to match the original implementation
def SinkSource(G, unknowns, max_iters=1000, delta=0.0001):
    """ Iterate over the graph until all of the scores of the unknowns have converged
        *max_iters*: maximum number of iterations
    """
    converged = False
    iters = 0
    start_time = time.time()

    for iters in trange(max_iters):
        for u in unknowns:
            data = G.node[u]
            # update the value of iterate all of the neighbors of u
            neighbor_sum = 0
            for v in G[u]:
                neighbor_sum += G[u][v]['weight'] * G.node[v]['prev_s']

            data['s'] = (neighbor_sum + data['c']) / data['degree']
            #data['s'] = neighbor_sum / data['degree']
            #G.node[u]['s'] = neighbor_sum / float(G.degree(u, weight='weight'))

        converged = True
        # Update the previous value, as well as how much the score changed by (delta)
        for u in unknowns:
            data = G.node[u]
            data['d'] = abs(data['prev_s'] - data['s'])
            # check if the score changed enough to indicate the algorithm hasn't yet converged
            if not converged or data['d'] >= delta:
                converged = False

            data['prev_s'] = data['s']

        #converged = hasConverged(G, unknowns, delta=delta)
        if converged:
            break

    total_time = time.time() - start_time
    print("SinkSource converged after %d iterations (%0.2f sec)" % (iters, total_time))
    #print("Scores of the top k nodes:")
    #print("\t%s" % (', '.join(["%s: %s" % (u, G.node[u]['s']) for u in list(unknowns)[:15]])))

    return G


# this is supposed to match the original implementation
def SinkSourceUB(G, unknowns, k=10, max_iters=1000, delta=0.0001):
    """
    Iterate over the graph until all of the scores of the unknowns have converged
        *max_iters*: maximum number of iterations
    *k*: number of nodes to get scores for
    """
    converged = False
    iters = 0
    start_time = time.time()

    for iters in trange(max_iters):
        for u in unknowns:
            data = G.node[u]
            # update the value of iterate all of the neighbors of u
            neighbor_sum = 0
            # upper bound (ub) sum is computed the same way
            ub_sum = 0
            for v in G[u]:
                w = G.edges[u,v]['weight']
                neighbor_sum += w * G.node[v]['prev_s']
                ub_sum += w * G.node[v]['prev_ub']

            data['s'] = (neighbor_sum + data['c']) / data['degree']
            data['ub'] = (ub_sum + data['c']) / data['degree']
            #G.node[u]['s'] = neighbor_sum / float(G.degree(u, weight='weight'))

        converged = True
        # Update the previous value, as well as how much the score changed by (delta)
        for u in unknowns:
            data = G.node[u]
            data['d'] = data['s'] - data['prev_s']
            data['d_ub'] = data['prev_ub'] - data['ub']
            # check if the score changed enough to indicate the algorithm hasn't yet converged
            if not converged or data['d'] >= delta:
                converged = False

            data['prev_s'] = data['s']
            data['prev_ub'] = data['ub']

        #converged = hasConverged(G, unknowns, delta=delta)
        if converged:
            break

        # see if there are any nodes to prune
        if iters > 10:
            # only the unknowns have the 's' and 'ub' attributes
            scores = {n:G.node[n]['s'] for n in unknowns}
            k_score = scores[sorted(scores, key=scores.get, reverse=True)[k]]
            ubs = {n:G.node[n]['ub'] for n in unknowns}

            nodes_to_prune = set()
            for u in sorted(ubs, key=ubs.get):
                if ubs[u] < k_score:
                    nodes_to_prune.add(u)
                else:
                    break

            if len(nodes_to_prune) > 0:
                tqdm.write("pruning %d nodes (out of %d) at iter: %d" % (len(nodes_to_prune), len(unknowns), iters))
                #pdb.set_trace()
                # remove them from the unknowns, and remove their current score
                unknowns.difference_update(nodes_to_prune)
                for u in nodes_to_prune:
                    data = G.node[u]
                    # fix their score
                    # TODO figure out the best way to fix the score
                    data['prev_s'] = data['s'] + (data['d'] / 2.0)
                    data['prev_ub'] = data['ub'] - (data['d_ub'] / 2.0)
                    #del data['s']
                    #del data['ub']
                    # Rather than fix this node's score, which will require its neighboring nodes to continue to perform the computations,
                    # Update the neighboring node's 'c' which is the amount of score recieved from fixed neighbors
                    # TODO this is not setup for directed graphs
                    for v in G[u]:
                        w = G.edges[u,v]['weight']
                        G.node[v]['c'] += w * data['prev_s']
                        G.node[v]['c_ub'] += w * data['prev_ub']
                G.remove_nodes_from(nodes_to_prune)

    total_time = time.time() - start_time
    print("SinkSource converged after %d iterations (%0.2f sec)" % (iters, total_time))
    #print("Scores of the top k nodes:")
    #print("\t%s" % (', '.join(["%s: %0.2f" % (u, G.node[u]['s']) for u in list(unknowns)[:k]])))

    return G


#def hasConverged(G, unknowns, delta=0.0001):
#    converged = True
#    for n in unknowns:
#        d = G.node[n]['d']
#        if d >= delta:
#            converged = False
#            break
#
#    return converged


def setupScores(G, positives, negatives, unknowns):
    """
    """

    print("Initializing scores and seteting up network")
    # initialize all of the scores of the nodes to 0
    for u in unknowns:
        data = G.node[u]
        data['s'] = 0
        data['prev_s'] = 0
        data['ub'] = 1
        data['prev_ub'] = 1
        data['degree'] = float(G.degree(u, weight='weight'))
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
        data['c'] = 0
        data['c_ub'] = 0
        for p in positives:
            if G.has_edge(u,p):
                data['c'] += G.edges[u,p]['weight']
                data['c_ub'] += G.edges[u,p]['weight']
        # don't add the negatives weight. Those only contribute to the denominator
        #for n in negatives:
        #    if G.has_edge(u,n):
        #        data['c'] += G.edges[u,n]['weight']

    # now remove the positive and negative nodes from the network
    G.remove_nodes_from(positives.union(negatives))

    return G


def setupSinkSourcePlus(G, unknowns, neg_node="n", neg_weight=1):
    """ Adds an extra negative node and connects it to all unknowns
    with a weight of *negative_weight*
    """

    G.add_weighted_edges_from([(n, neg_node, neg_weight) for n in unknowns])

    return G


def runSinkSource(G, positives, negatives=None, k=None, max_iters=1000, delta=0.0001):
    """
    *positives*: set of positive nodes
    *negeatives*: set of negative nodes
    *k*: Run with upper and lower bounds to prune unnecessary nodes
    """
    sinksourceplus = False
    if negatives is None:
        sinksourceplus = True
        negatives = set(['n'])
    unknowns = set(G.nodes()).difference(positives).difference(negatives)

    if sinksourceplus:
        # set the negative node's weight to be the average of the weights connected to the positive nodes
        neg_weight = sum([G.degree(p, weight='weight')/G.degree(p) for p in positives]) / float(len(positives))
        print("Setting the weight of the negative node to the avg of the positive nodes: %s" % (str(neg_weight)))
        G = setupSinkSourcePlus(G, unknowns, neg_node='n', neg_weight=neg_weight)

    setupScores(G, positives, negatives, unknowns)

    # replace the nodes with integers to speed things up
    index = 1
    node2int = {}
    int2node = {}
    for n in G.nodes():
        node2int[n] = index
        int2node[index] = n
        index += 1
    # see also convert_node_labels_to_integers
    G = nx.relabel_nodes(G,node2int)
    # the positives and negatives were already removed from the graph,
    # so everything left is an unknown
    unknowns = set(G.nodes())

    if k is not None:
        print("\tGetting top %d predictions" % (k))
        G = SinkSourceUB(G, unknowns, k=k, max_iters=max_iters, delta=delta)
    else:
        G = SinkSource(G, unknowns, max_iters=max_iters, delta=delta)

    # switch the nodes back to their names
    G = nx.relabel_nodes(G,int2node)

    return G
