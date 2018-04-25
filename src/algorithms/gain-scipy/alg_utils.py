
from scipy import sparse
from scipy.sparse import csr_matrix
import numpy as np
import networkx as nx


def convert_nodes_to_int(G):
    index = 0
    node2int = {}
    int2node = {}
    for n in G.nodes():
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


def setupScores(P, positives, negatives=None, a=1):
    """
    """
    #print("Initializing scores and setting up network")
    pos_vec = np.zeros(P.shape[0])
    pos_vec[positives] = 1
    #if negatives is not None:
    #    pos_vec[negatives] = -1

    # f contains the fixed amount of score coming from positive nodes
    f = a*csr_matrix.dot(P, pos_vec)

    if negatives is None:
        fixed_nodes = positives
    else:
        fixed_nodes = np.concatenate([positives, negatives])

    # keep track of the original node integers 
    # to be able to map back to the original node names
    int2int = {}
    index_diff = 0
    for i in range(len(f)):
        if i in fixed_nodes:
            index_diff += 1
            continue
        int2int[i - index_diff] = i

    # remove the fixed nodes from the graph
    # removing the fixed nodes is slightly faster than selecting the unknown rows
    P = delete_nodes(P, fixed_nodes)
    # and from f
    f = np.delete(f, fixed_nodes)
    assert P.shape[0] == P.shape[1], "Matrix is not square"
    assert P.shape[1] == len(f), "f doesn't match size of P"

    return P, f, int2int


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
