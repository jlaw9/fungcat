
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


def nonReachableNodes(G, nodes):
    """
    Returns a set of nodes that cannot be reached from the input set of (positive) nodes
    """
    curr_nodes = nodes
    curr_neighbors = set()
    not_visited = set(G.nodes())
    # BFS to every node reachable from a positive
    while len(curr_nodes) > 0:
        # continue iterating until there are no more nodes
        for u in curr_nodes:
            curr_neighbors.update(set(G.neighbors(u)) & not_visited)
        not_visited.difference_update(curr_nodes)
        curr_nodes = curr_neighbors
        curr_neighbors = set()
    return not_visited
