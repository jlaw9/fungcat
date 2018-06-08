
import os, sys
from scipy import sparse
from scipy.sparse import csr_matrix, csgraph
import numpy as np
import networkx as nx
from collections import defaultdict
import time
sys.path.append("src")
import utils.file_utils as utils


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


#def normalizeGraphEdgeWeights(W):
def normalizeGraphEdgeWeights(W, l=None):
    """
    *W*: weighted network as a scipy sparse matrix in csr format
    *l*: SinkSourcePlus lambda parameter
    """
    # normalize the matrix
    # by dividing every edge by the node's degree (row sum)
    if l is None:
        P = W.multiply(csr_matrix(1/W.sum(axis=1).astype(float)))
    else:
        P = W.multiply(csr_matrix(1/(l+W.sum(axis=1).astype(float))))
    return P
    # this gives a memory error likely because it first converts the matrix to a numpy matrix
    #return W / W.sum(axis=1)
#def normalizeGraphEdgeWeights(W):
#    """
#    *W*: weighted network as a scipy sparse matrix in csr format
#    """
#    # normalizing the matrix
#    # sum the weights in the row to get the degree of each node
#    deg = np.asarray(W.sum(axis=1)).flatten()
#    deg = np.divide(1., np.sqrt(deg))
#    deg[np.isinf(deg)] = 0
#    D = sparse.diags(deg).tocsr()
#    # normalize W by multiplying D^(-1/2) * W * D^(-1/2)
#    W = csr_matrix.dot(D, csr_matrix.dot(W, D))
#    return W


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



def check_fixed_rankings(LBs, UBs, unranked_nodes, nodes_to_rank=None):
    """
    *nodes_to_rank*: a set of nodes for which to check which nodes have an overlapping UB/LB.
        In other words, get the nodes that are causing the given set of nodes to not have their ranks fixed
    """
    # Create two lists to sort.
    # One with the string "LB" and one with the string "UB".
    # Combine both lists. If for each node the 'LB' and 'UB' strings are next to each other, then the nodes rank is fixed
    # TODO use the epsUB parameter
    # find all of the nodes in the top-k whose rankings are fixed
    # nlogn comparisons
    all_scores = [(n, 'LB', LBs[n]) for n in unranked_nodes]
    all_scores += [(n, 'UB', UBs[n]) for n in unranked_nodes]
    # sort first by the score, and then by the node name which will break ties
    # TODO this is really slow at the start. 
    # Maybe I could speed it up by starting with a random sampling of nodes, 
    # and once those are ranked move up to all nodes
    all_scores_sorted = sorted(all_scores, key=lambda x: (x[2], x[0]), reverse=True)
    i = 0
    fixed_nodes = set()
    if nodes_to_rank is None:
        while i+1 < len(all_scores_sorted):
            if all_scores_sorted[i][0] == all_scores_sorted[i+1][0]:
                # temp check
                if all_scores_sorted[i][2] == 'LB':
                    print("Warning: %s LB < UB: %s < %s" % (all_scores_sorted[i][0],
                        all_scores_sorted[i][2], all_scores_sorted[i+1][2])) 
                fixed_nodes.add(all_scores_sorted[i][0])
                i += 1
            i += 1
    else:
        conflicting_nodes = set()
        conflicted_node = -1
        conflicts = False
        while i < len(all_scores_sorted):
            curr_n = all_scores_sorted[i][0]
            if i+1 < len(all_scores_sorted) and curr_n == all_scores_sorted[i+1][0]:
                fixed_nodes.add(curr_n)
                i += 1
            # if nodes_to_rank is specified,
            # then check to see which nodes are conflicting with the ranks
            elif curr_n in nodes_to_rank and all_scores_sorted[i][1] == "UB":
                conflicts = True
                conflicted_node = curr_n 
                conflicting_nodes.add(curr_n) 
            # if we've hit the LB, then stop checking
            elif curr_n == conflicted_node:
                conflicts = False
            elif conflicts is True:
                conflicting_nodes.add(curr_n)
            i += 1

        #unranked_nodes = conflicting_nodes
        return conflicting_nodes

    # n^2 comparisons
#                fixed_nodes = unranked_nodes.copy()
#                for u in unranked_nodes:
#                    for v in unranked_nodes:
#                        # don't check the same pair of nodes twice
#                        if v <= u:
#                            continue
#                        # if the LB and UB of these two nodes overlap, then these two node's rankings are not fixed
#                        if LBs[u] < UBs[v] - self.epsUB and \
#                                UBs[u] - self.epsUB > LBs[v]:
#                            fixed_nodes.discard(u)
#                            fixed_nodes.discard(v)
##                            if len(unranked_nodes) < 70:
##                                print("\t\tLBs[u] < UBs[v] and UBs[u] > LBs[v]: %0.6f < %0.6f, %0.6f > %0.6f" % (LBs[u], UBs[v], UBs[u], LBs[v]))
#                            break

    return fixed_nodes


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

## copied from here: https://stackoverflow.com/a/26504995
#def delete_rows_csr(mat, indices):
#    """
#    Remove the rows denoted by ``indices`` form the CSR sparse matrix ``mat``.
#    """
#    if not isinstance(mat, scipy.sparse.csr_matrix):
#        raise ValueError("works only for CSR format -- use .tocsr() first")
#    indices = list(indices)
#    mask = np.ones(mat.shape[0], dtype=bool)
#    mask[indices] = False
#    return mat[mask]
#
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

