
import networkx as nx


def convert_labels_to_int(G):
    index = 0
    node2int = {}
    int2node = {}
    for n in G.nodes():
        node2int[n] = index
        int2node[index] = n
        index += 1
    # see also convert_node_labels_to_integers
    H = nx.relabel_nodes(G,node2int)
    return H, node2int, int2node


def normalizeGraphEdgeWeights(G):
    """
    Row-normalize the graph G.
    Returns a directed G (shallow copy)
    """
    H = nx.DiGraph(G)
    #H = G.to_directed()
    for u in H.nodes():
        deg = float(H.out_degree(u, weight='weight'))
        for v, data in H.adj[u].items():
            data['weight'] = data['weight'] / deg

    return H


