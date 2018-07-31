
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
    H = nx.relabel_nodes(G,node2int, copy=False)
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


# TODO implement a better method for getting the goid annotations
#def parse_pos_neg_matrix(pos_neg_file, include_negatives=True):
#
#    # TODO figure out a better way to parse these files
#    goid_prots_file = pos_neg_file.replace(".tsv", "-pos-prots.txt")
#    goid_neg_file = pos_neg_file.replace(".tsv", "-neg-prots.txt")
#    all_goid_prots = {}
#    all_goid_neg = {}
#    #if forced is False and os.path.isfile(goid_prots_file):
#    if os.path.isfile(goid_prots_file) and (include_negatives is False or os.path.isfile(goid_neg_file)):
#        print("Reading goid protein annotations from %s" % (goid_prots_file))
#        all_goid_prots = {goid: set(prots.split(',')) for goid, prots in utils.readColumns(goid_prots_file, 1, 2)}
#        if include_negatives is True:
#            all_goid_neg = {goid: set(prots.split(',')) for goid, prots in utils.readColumns(goid_neg_file, 1, 2)}
#
#    else:
#        pos_neg_files = ['%s/nonieapos-neg-bp-%d.tsv' % (in_dir, cutoff), '%s/nonieapos-neg-mf-%d.tsv' % (in_dir, cutoff)]
#        for pos_neg_file in pos_neg_files:
#            goid_prots, goid_neg = parse_pos_neg_matrix(pos_neg_file)
#            all_goid_prots.update(goid_prots)
#            all_goid_neg.update(goid_neg)
#
#        with open(goid_prots_file, 'w') as out:
#            out.write(''.join(["%s\t%s\n" % (goid, ','.join(prots)) for goid, prots in all_goid_prots.items()]))
#
#        with open(goid_neg_file, 'w') as out:
#            out.write(''.join(["%s\t%s\n" % (goid, ','.join(prots)) for goid, prots in all_goid_neg.items()]))
#    print("Reading positive and negative annotations for each protein from %s" % (pos_neg_file))
#    goid_prots = defaultdict(set)
#    goid_neg = defaultdict(set)
#    num_prots = 0
#    #pos_neg = df.read_csv(pos_neg_file)
#    # for each GO term, get the set of positives and negatives
#    with open(pos_neg_file, 'r') as f:
#        goterms = f.readline().rstrip().split('\t')[1:]
#        for line in f:
#            line = line.rstrip().split('\t')
#            num_prots += 1
#            prot = line[0]
#            vals = line[1:]
#            for i in range(0, len(vals)):
#                val = int(vals[i])
#                if val == 1:
#                    goid_prots[goterms[i]].add(prot)
#                elif val == -1:
#                    goid_neg[goterms[i]].add(prot)
#
#    print("\t%d GO terms, %d prots" % (len(goid_prots), num_prots))
#
#    return goid_prots, goid_neg

