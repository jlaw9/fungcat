
import os, sys
from scipy import sparse
from scipy.sparse import csr_matrix, csgraph
import numpy as np
import networkx as nx
from collections import defaultdict
import time
from tqdm import tqdm
sys.path.append("src")
import utils.file_utils as utils
# needed for evaluation metrics
#try:
from sklearn import metrics
#except ImportError:
#    pass


def select_goterms(only_functions_file=None, goterms=None):
    selected_goterms = set()
    if only_functions_file is not None:
        selected_goterms = utils.readItemSet(only_functions_file, 1)
        # if the leading GO ID isn't present, then add it
        if list(selected_goterms)[0][0:3] != "GO:":
            selected_goterms = set(["GO:" + "0"*(7-len(str(x))) + str(x) for x in selected_goterms])
    goterms = set() if goterms is None else set(goterms)
    selected_goterms.update(goterms)
    if len(selected_goterms) == 0:
        selected_goterms = None
    return selected_goterms


def parse_pos_neg_files(pos_neg_files, goterms=None):
    # get the positives and negatives from the matrix
    all_goid_prots = {}
    all_goid_neg = {}
    if len(pos_neg_files) == 1 and pos_neg_files[0] == '-':
        print("Using GAIN annotations instead of pos_neg_file")
        # TODO compare the predictions made by GAIN and my implementation
        all_goid_prots = parse_gain_file(f_settings.GOA_ALL_FUN_FILE_NOIEA)
        all_goid_neg = {goid: set() for goid in all_goid_prots} 
    else:
        for pos_neg_file in pos_neg_files:
            #goid_prots, goid_neg = self.parse_pos_neg_matrix(self.pos_neg_file)
            goid_prots, goid_neg = parse_pos_neg_file(pos_neg_file, goterms=goterms)
            all_goid_prots.update(goid_prots)
            all_goid_neg.update(goid_neg)

    return all_goid_prots, all_goid_neg


def parse_pos_neg_file(pos_neg_file, goterms=None):
    print("Reading positive and negative annotations for each protein from %s" % (pos_neg_file))
    goid_prots = {}
    goid_neg = {}
    all_prots = set()
    # TODO possibly use pickle
    if not os.path.isfile(pos_neg_file):
        print("Warning: %s file not found" % (pos_neg_file))
        return goid_prots, goid_neg

        #for goid, pos_neg_assignment, prots in utils.readColumns(pos_neg_file, 1,2,3):
    with open(pos_neg_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            goid, pos_neg_assignment, prots = line.rstrip().split('\t')[:3]
            if goterms and goid not in goterms:
                continue
            prots = set(prots.split(','))
            if int(pos_neg_assignment) == 1:
                goid_prots[goid] = prots
            elif int(pos_neg_assignment) == -1:
                goid_neg[goid] = prots

            all_prots.update(prots)

    print("\t%d GO terms, %d prots" % (len(goid_prots), len(all_prots)))

    return goid_prots, goid_neg


def parse_gain_file(gain_file, goterms=None):
    print("Reading annotations from GAIN file %s. Assuming annotations have already been propogated up the GO DAG" % (gain_file))
    goid_prots = defaultdict(set)

    with open(gain_file, 'r') as f:
        # columns: 'orf', 'goid', 'hierarchy', 'evidencecode', 'annotation type' (described here: http://bioinformatics.cs.vt.edu/~murali/software/biorithm/gain.html)
        header_line = f.readline().rstrip().split('\t')
        # *orf* A systematic name for the gene.
        orf_index = header_line.index('orf')
        # *goid* The ID of the GO function. You can leave in the "GO:0+" prefix for a function. GAIN will strip it out.
        goid_index = header_line.index('goid')
        goname_index = header_line.index('goname')
        # *hierarchy* The GO category the function belongs to ("c", "f", and "p")
        hierarchy_index = header_line.index('hierarchy')
        # *evidencecode* The evidence code for an annotation
        evidencecode_index = header_line.index('evidencecode')
        # *annotation type* The value is either 1 (annotated), 0 (unknown), or -1 (not annotated)
        ann_index = header_line.index('annotation type')

        for line in f:
            line = line.rstrip().split('\t')
            prot = line[orf_index]
            goid = line[goid_index]
            # convert the goterm id to the full ID
            goid = "GO:" + "0"*(7-len(goid)) + goid
            goname = line[goname_index]
            hierarchy = line[hierarchy_index]
            evidencecode = line[evidencecode_index]
            annotation_type = line[ann_index]

            # only track the annotations that are positive
            if annotation_type != "1":
                continue

            #hierarchies[goid] = hierarchy
            #gonames[goid] = goname
            if goterms is None or goid in goterms:
                goid_prots[goid].add(prot)

    return goid_prots


def convert_nodes_to_int(G):
    index = 0
    node2int = {}
    int2node = {}
    for n in sorted(G.nodes()):
        node2int[n] = index
        int2node[index] = n
        index += 1
    # see also convert_node_labels_to_integers
    G = nx.relabel_nodes(G,node2int, copy=False)
    return G, node2int, int2node
#def convert_nodes_to_int(nodes):
#    index = 0
#    node2int = {}
#    int2node = {}
#    for n in nodes:
#        node2int[n] = index
#        int2node[index] = n
#        index += 1
#    # see also convert_node_labels_to_integers
#    #H = nx.relabel_nodes(G,node2int, copy=False)
#    return node2int, int2node


def setup_sparse_network(network_file, node2idx_file=None, forced=False):
    """
    Takes a network file and converts it to a sparse matrix
    """
    sparse_net_file = network_file.replace('.txt', '.npz')
    if node2idx_file is None:
        node2idx_file = sparse_net_file + "-node-ids.txt"
    if forced is False and (os.path.isfile(sparse_net_file) and os.path.isfile(node2idx_file)):
        print("Reading network from %s" % (sparse_net_file))
        W = sparse.load_npz(sparse_net_file)
        print("\t%d nodes and %d edges" % (W.shape[0], len(W.data)/2))
        print("Reading node names from %s" % (node2idx_file))
        node2idx = {n: int(n2) for n, n2 in utils.readColumns(node2idx_file, 1, 2)}
        idx2node = {n2: n for n, n2 in node2idx.items()}
        prots = [idx2node[n] for n in sorted(idx2node)]
    elif os.path.isfile(network_file):
        print("Reading network from %s" % (network_file))
#        G = nx.Graph()
#        with open(network_file, 'r') as f:
#            for line in f:
#                if line[0] == "#":
#                    continue
#                u,v,w = line.rstrip().split('\t')[:3]
#                G.add_edge(u,v,weight=float(w))
#                #edges[(u,v)] = float(w)
#                #nodes.update(set([u,v]))
#        print("\t%d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))
#
#        print("\trelabeling node IDs with integers")
#        G, node2idx, idx2node = convert_nodes_to_int(G)
#        print("\twriting node2idx labels to %s" % (node2idx_file))
#        nodes_ids = sorted(idx2node)
#        with open(node2idx_file, 'w') as out:
#            out.write(''.join(["%s\t%s\n" % (idx2node[n], n) for n in nodes_ids]))
#
#        print("\tconverting to a scipy sparse matrix")
#        W = nx.to_scipy_sparse_matrix(G, nodelist=nodes_ids)
        # load the first three columns into vectors
        #u, v, val = np.loadtxt(filename).T[:3]
        u,v,w = [], [], []
        # TODO make sure the network is symmetrical
        with open(network_file, 'r') as f:
            # add tqdm?
            for line in tqdm(f, total=120000000):
                line = line.rstrip().split('\t')
                u.append(line[0])
                v.append(line[1])
                w.append(float(line[2]))
        print("\tconverting uniprot ids to node indexes / ids")
        # first convert the uniprot ids to node indexes / ids
        prots = sorted(set(list(u)) | set(list(v)))
        node2idx = {prot: i for i, prot in enumerate(prots)}
        i = [node2idx[n] for n in u]
        j = [node2idx[n] for n in v]
        print("\tcreating sparse matrix")
        #print(i,j,w)
        W = sparse.coo_matrix((w, (i, j)), shape=(len(prots), len(prots))).tocsr()
        # make sure it is symmetric
        if (W.T != W).nnz == 0:
            pass
        else:
            print("### Matrix not symmetric!")
            W = W + W.T
            print("### Matrix converted to symmetric.")
        #name = os.path.basename(net_file)
        print("\twriting sparse matrix to %s" % (sparse_net_file))
        sparse.save_npz(sparse_net_file, W)
        print("\twriting node2idx labels to %s" % (node2idx_file))
        with open(node2idx_file, 'w') as out:
            out.write(''.join(["%s\t%d\n" % (prot,i) for i, prot in enumerate(prots)]))
    else:
        print("Network %s not found. Quitting" % (network_file))
        sys.exit(1)

    return W, prots


def normalizeGraphEdgeWeights(W, ss_lambda=None, axis=1):
    """
    *W*: weighted network as a scipy sparse matrix in csr format
    *ss_lambda*: SinkSourcePlus lambda parameter
    *axis*: The axis to normalize by. 0 is columns, 1 is rows
    """
    # normalize the matrix
    # by dividing every edge by the node's degree (row sum)
    deg = np.asarray(W.sum(axis=axis)).flatten()
    if ss_lambda is None:
        deg = np.divide(1., deg)
    else:
        deg = np.divide(1., ss_lambda + deg)
    deg[np.isinf(deg)] = 0
    P = W.multiply(deg)
    return P.asformat(W.getformat())
    # this gives a memory error likely because it first converts the matrix to a numpy matrix
    #return W / W.sum(axis=1)


def _net_normalize(W, axis=0):
    """
    Normalize W by multiplying D^(-1/2) * W * D^(-1/2)
    This is used for GeneMANIA
    *W*: weighted network as a scipy sparse matrix in csr format
    """
    # normalizing the matrix
    # sum the weights in the columns to get the degree of each node
    deg = np.asarray(W.sum(axis=axis)).flatten()
    deg = np.divide(1., np.sqrt(deg))
    deg[np.isinf(deg)] = 0
    D = sparse.diags(deg)
    # normalize W by multiplying D^(-1/2) * W * D^(-1/2)
    P = D.dot(W.dot(D))
    return P


def get_goid_pos_neg(ann_matrix, i):
    """
    The matrix should be lil format as others don't have the getrowview option
    """
    # get the row corresponding to the current goids annotations 
    #goid_ann = ann_matrix[i,:].toarray().flatten()
    #positives = np.where(goid_ann > 0)[0]
    #negatives = np.where(goid_ann < 0)[0]
    # may be faster with a lil matrix, but takes a lot more RAM
    #goid_ann = ann_matrix.getrowview(i)
    goid_ann = ann_matrix[i,:]
    positives = (goid_ann > 0).nonzero()[1]
    negatives = (goid_ann < 0).nonzero()[1]
    return positives, negatives


def setupScores(P, positives, negatives=None, a=1, remove_nonreachable=True, verbose=False):
    """
    """
    #print("Initializing scores and setting up network")
    pos_vec = np.zeros(P.shape[0])
    pos_vec[positives] = 1
    #if negatives is not None:
    #    pos_vec[negatives] = -1

    # f contains the fixed amount of score coming from positive nodes
    f = a*csr_matrix.dot(P, pos_vec)

    if remove_nonreachable is True:
        node2idx, idx2node = {}, {}
        # remove the negatives first and then remove the non-reachable nodes
        if negatives is not None:
            node2idx, idx2node = build_index_map(range(len(f)), negatives)
            P = delete_nodes(P, negatives)
            f = np.delete(f, negatives)
            #fixed_nodes = np.concatenate([positives, negatives])
            positives = set(node2idx[n] for n in positives)
        positives = set(list(positives))
        fixed_nodes = positives 

        start = time.time()
        # also remove nodes that aren't reachable from a positive 
        # find the connected_components. If a component doesn't have a positive, then remove the nodes of that component
        num_ccs, node_comp = csgraph.connected_components(P, directed=False)
        # build a dictionary of nodes in each component
        ccs = defaultdict(set)
        # check to see which components have a positive node in them
        pos_comp = set()
        for n in range(len(node_comp)):
            comp = node_comp[n]
            ccs[comp].add(n)
            if comp in pos_comp:
                continue
            if n in positives:
                pos_comp.add(comp)

        non_reachable_ccs = set(ccs.keys()) - pos_comp
        not_reachable_from_pos = set(n for cc in non_reachable_ccs for n in ccs[cc])
#        # use matrix multiplication instead
#        reachable_nodes = get_reachable_nodes(P, positives)
#        print(len(reachable_nodes), P.shape[0] - len(reachable_nodes))
        if verbose:
            print("%d nodes not reachable from a positive. Removing them from the graph" % (len(not_reachable_from_pos)))
            print("\ttook %0.4f sec" % (time.time() - start))
        # combine them to be removed
        fixed_nodes = positives | not_reachable_from_pos

        node2idx2, idx2node2 = build_index_map(range(len(f)), fixed_nodes)
        if negatives is not None:
            # change the mapping to be from the deleted nodes to the original node ids
            node2idx = {n: node2idx2[node2idx[n]] for n in node2idx if node2idx[n] in node2idx2}
            idx2node = {node2idx[n]: n for n in node2idx}
        else:
            node2idx, idx2node = node2idx2, idx2node2 
    else:
        fixed_nodes = positives 
        if negatives is not None:
            fixed_nodes = np.concatenate([positives, negatives])
        node2idx, idx2node = build_index_map(range(len(f)), set(list(fixed_nodes)))
    # removing the fixed nodes is slightly faster than selecting the unknown rows
    # remove the fixed nodes from the graph
    fixed_nodes = np.asarray(list(fixed_nodes)) if not isinstance(fixed_nodes, np.ndarray) else fixed_nodes
    P = delete_nodes(P, fixed_nodes)
    # and from f
    f = np.delete(f, fixed_nodes)
    assert P.shape[0] == P.shape[1], "Matrix is not square"
    assert P.shape[1] == len(f), "f doesn't match size of P"

    return P, f, node2idx, idx2node


def build_index_map(nodes, nodes_to_remove):
    """
    *returns*: a dictionary of the original node ids/indices to the current indices, as well as the reverse
    """
    # keep track of the original node integers 
    # to be able to map back to the original node names
    node2idx = {}
    idx2node = {}
    index_diff = 0
    for i in nodes:
        if i in nodes_to_remove:
            index_diff += 1
            continue
        idx2node[i - index_diff] = i
        node2idx[i] = i - index_diff

    return node2idx, idx2node 


#def get_reachable_nodes(P, start_nodes):
#    """
#    Get the nodes that are reachable from a set of starting nodes
#    *returns*: a np array of reachable nodes
#    """
#    start = np.zeros(P.shape[0])
#    start[list(start_nodes)] += 1
#    #Pbool = P.astype(bool)
#    #visited = sparse.diags(start).astype(bool)
#    visited = sparse.diags(start)
#    # R starts out as a matrix of zeros
#    R = sparse.csr_matrix(P.shape)
#    while visited.count_nonzero() > 0:
#        R = R + visited
#        # TODO the graph is undirected so it keeps repeating cycles
#        # I need a way to not follow cycles
#        visited = visited.dot(P)
#        print(visited.count_nonzero())
#    # return all of the reachable nodes
#    reachable_nodes = set(R.getnnz(axis=0).nonzero()[0])
#    return reachable_nodes


def check_fixed_rankings(LBs, UBs, unranked_nodes=None, unr_pos_nodes=None, unr_neg_nodes=None):
    """
    *nodes_to_rank*: a set of nodes for which to check which nodes have an overlapping UB/LB.
        In other words, get the nodes that are causing the given set of nodes to not have their ranks fixed
    UPDATE: 
    *unr_pos_nodes*: set of positive nodes that are not fixed. 
        only need to check for overlap with negative nodes 
    *unr_neg_nodes*: set of negative nodes that are not fixed. 
        only need to check for overlap with positive nodes
    """
    # TODO use the epsUB parameter
    # find all of the nodes in the top-k whose rankings are fixed
    # n comparisons
    all_scores = []
    # also keep track of the # of nodes in a row that have overlapping upper or lower bounds
    # for now just keep track of the biggest
    max_unranked_stretch = 0
    i = 0
    if unranked_nodes is not None:
        for n in unranked_nodes:
            all_scores.append((n, LBs[n]))
        all_scores_sorted = sorted(all_scores, key=lambda x: (x[1]), reverse=True)
        # the fixed nodes are the nodes that are not in the 
        # "still not fixed" set
        still_not_fixed_nodes = set()
        # for every node, check if the next node's LB+UB > the curr node's LB.
        # If so, the node is not yet fixed
        while i+1 < len(all_scores_sorted):
            curr_LB = all_scores_sorted[i][1]
            curr_i = i
            while i+1 < len(all_scores_sorted) and \
                curr_LB < UBs[all_scores_sorted[i+1][0]]:
                still_not_fixed_nodes.add(all_scores_sorted[i+1][0])
                #print("i+1: %d not fixed" % (i+1))
                i += 1
            if curr_i != i:
                #print("i: %d not fixed" % (curr_i))
                still_not_fixed_nodes.add(all_scores_sorted[curr_i][0])
                if i - curr_i > max_unranked_stretch:
                    max_unranked_stretch = i - curr_i
            if curr_i == i:
                i += 1
            #    fixed_nodes.add(all_scores_sorted[i][0])
        return still_not_fixed_nodes, max_unranked_stretch
    elif unr_pos_nodes is not None and unr_neg_nodes is not None:
        # if there aren't any nodes to check their ranking, then simply return
        if len(unr_pos_nodes) == 0 or len(unr_neg_nodes) == 0:
            return set(), set()
        for n in unr_pos_nodes:
            all_scores.append((n, LBs[n]))
        for n in unr_neg_nodes:
            all_scores.append((n, LBs[n]))
        all_scores_sorted = sorted(all_scores, key=lambda x: (x[1]), reverse=True)
        fixed_nodes = unr_pos_nodes | unr_neg_nodes
        # for every node, check if the next node's LB+UB > the curr node's LB.
        # and if one of the overlapping nodes is opposite 
        # If so, the node is not yet fixed
        curr_node = all_scores_sorted[0][0]
        opp_set_pos = False if curr_node in unr_pos_nodes else True
        opp_set = unr_pos_nodes if opp_set_pos else unr_neg_nodes
        while i+1 < len(all_scores_sorted):
            curr_node = all_scores_sorted[i][0]
            curr_LB = all_scores_sorted[i][1]
            # if this is a positive, just check the negatives
            # and vice versa
            curr_i = i
            opp_overlap = False 
            last_opp_node = None
            while i+1 < len(all_scores_sorted) and \
                curr_LB < UBs[all_scores_sorted[i+1][0]]:
                next_node = all_scores_sorted[i+1][0]
                if next_node in opp_set:
                    opp_overlap = True
                    last_opp_node = i+1 
                i += 1
            if opp_overlap is True:
                # if there was an overlap with an opposite node,
                # all of these are not fixed
                for j in range(curr_i, i+1):
                    j_node = all_scores_sorted[j][0]
                    fixed_nodes.discard(j_node)
                i = last_opp_node
                # flip the opposite set
                opp_set_pos = False if opp_set_pos else True
                opp_set = unr_pos_nodes if opp_set_pos else unr_neg_nodes
            else:
                # only need to increment here if there was no overlap
                i += 1
                
        unr_pos_nodes -= fixed_nodes 
        unr_neg_nodes -= fixed_nodes 
        return unr_pos_nodes, unr_neg_nodes
    else:
        print("Error: need to pass either the 'unranked_nodes' set or both 'pos_nodes' and 'neg_nodes'")
        return


# this is like 10x slower than python's sort function
def insertionSort(alist, tuple_idx=None):
    if tuple_idx is None:
        for index in range(1,len(alist)):
            currentvalue = alist[index]
            position = index

            while position > 0 and alist[position-1] > currentvalue:
                alist[position] = alist[position-1]
                position = position-1

            alist[position] = currentvalue
    else:
        for index in range(1,len(alist)):
            currentvalue = alist[index]
            position = index

            while position > 0 and alist[position-1][tuple_idx] > currentvalue[tuple_idx]:
                alist[position] = alist[position-1]
                position = position-1

            alist[position] = currentvalue


def get_neighbors(P):
    neighbors = [[] for k in range(P.shape[0])]
    rows, cols = P.nonzero()
    for i in range(len(rows)):
        neighbors[rows[i]].append(cols[i])
    #neighbors = [set(x) for x in neighbors]
    neighbors = np.array([set(x) for x in neighbors])
    return neighbors


def neighbors(P, n):
    """ return the indices in the row of n (outgoing neighbors)
    that have a non-zero weight
    """
    # tooo sloww
    #return P[n].nonzero()[1]
    return P.indices[P.indptr[n]:P.indptr[n+1]]


def delete_nodes(mat, indices):
    """
    Remove the rows and columns denoted by ``indices`` form the CSR sparse matrix ``mat``.
    """
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask, :][:, mask]


# slightly faster than mat[indices, :][:, indices]
def select_nodes(mat, indices):
    """
    Select the rows and columns denoted by ``indices`` form the CSR sparse matrix ``mat``.
    Equivalent to getting a subnetwork of a graph
    """
    mask = np.zeros(mat.shape[0], dtype=bool)
    mask[indices] = True
    return mat[mask, :][:, mask]


# copied from here: https://stackoverflow.com/a/26504995
def delete_rows_csr(mat, indices):
    """
    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
    """
    if not isinstance(mat, sparse.csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")
    indices = list(indices)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[indices] = False
    return mat[mask]

#def delete_cols_csc(mat, indices):
#    """
#    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
#    """
#    assert isinstance(mat, scipy.sparse.csc_matrix), \
#            "works only for CSR format -- use .tocsr() first" 
#    indices = list(indices)
#    mask = np.ones(mat.shape[1], dtype=bool)
#    #print(mask)
#    #print(len(mask))
#    mask[indices] = False
#    return mat[:,mask]


def compute_eval_measures(scores, positives, negatives=None, track_pos_neg=False):
    """
    Compute the precision and false-positive rate at each change in recall (true-positive rate)
    *scores*: dictionary containing a score for each node
    *negatives*: if negatives are given, then the FP will only be from the set of negatives given
    *track_pos_neg*: if specified, track the score and rank of the positive and negative nodes,
        and return a tuple of the node ids in order of their score, their score, their idx, and 1/-1 for pos/neg
    """
    #f1_score = metrics.f1score(positives, 
    #num_unknowns = len(scores) - len(positives) 
    positives = set(positives)
    check_negatives = False
    if negatives is not None:
        check_negatives = True 
        negatives = set(negatives)
    else:
        print("TODO. Treating all non-positives as negatives not yet implemented.")
    # compute the precision and recall at each change in recall
    # TODO I should call numpy argsort to ensure I'm using the full precision when comparing values
    # use np.argsort
    #nodes_sorted_by_scores = sorted(scores, key=scores.get, reverse=True)
    nodes_sorted_by_scores = np.argsort(scores)[::-1]
    #print("computing the rank of positive nodes")
    # this is really slow...
    #pos_ranks = sorted([nodes_sorted_by_scores.index(p)+1 for p in positives])
    #print("%d positives, %d pos_ranks" % (len(positives), len(pos_ranks)))
    #print(pos_ranks)
    #print([scores[s] for s in nodes_sorted_by_scores[:pos_ranks[0]+1]])
    precision = [1]
    recall = [0]
    fpr = []
    pos_neg_stats = []  # tuple containing the node, score and idx
    # TP is the # of correctly predicted positives so far
    TP = 0
    FP = 0
    rec = 0
    for i, n in enumerate(nodes_sorted_by_scores):
        # TODO this could be slow if there are many positives
        if n in positives:
            TP += 1
            # precisions is the # of true positives / # true positives + # of false positives (or the total # of predictions)
            precision.append(TP / float(TP + FP))
            # recall is the # of recovered positives TP / TP + FN (total # of positives)
            rec = TP / float(len(positives))
            recall.append(rec)
            # fpr is the FP / FP + TN
            fpr.append((rec, FP / float(len(negatives))))
            if track_pos_neg:
                pos_neg_stats.append((n, scores[n], i, 1)) 
        elif check_negatives is False or n in negatives:
            FP += 1
            fpr.append((rec, FP / float(len(negatives))))
        #else:
        #    continue

    # TODO how should I handle this case?
    if len(precision) == 0:
        precision.append(0)
        recall.append(1)

    #print(precision[0], recall[0], fpr[0])

    if track_pos_neg:
        return precision, recall, fpr, pos_neg_stats
    else:
        return precision, recall, fpr




def compute_fmax(prec, rec):
    f_measures = []
    for i in range(len(prec)):
        p, r = prec[i], rec[i]
        if p+r == 0:
            harmonic_mean = 0
        else:
            harmonic_mean = (2*p*r)/(p+r)
        f_measures.append(harmonic_mean)
    return max(f_measures)

def compute_avgp(prec, rec):
    # average precision score
    # see http://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html#sklearn.metrics.average_precision_score
    avgp = 0
    prev_r = 0 
    for p,r in zip(prec, rec):
        recall_change = r - prev_r
        avgp += (recall_change*p)
        prev_r = r
    #avgp = avgp / float(len(alg_prec_rec))
    return avgp

def compute_auprc(prec, rec):
    auprc = metrics.auc(rec, prec)
    return auprc

def compute_auroc(tpr, fpr):
    auroc = metrics.auc(fpr, tpr)
    return auroc
