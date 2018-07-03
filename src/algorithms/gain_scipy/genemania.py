
# Python implementation of SinkSource
import time
#from tqdm import trange, tqdm
import alg_utils
import numpy as np
#mport gc
from scipy.sparse import csr_matrix, eye, diags
from scipy.sparse.linalg import spsolve, cg


def setup_laplacian(W):
    # first normalize the network
    P = alg_utils._net_normalize(W)
    # take the column sum and set them as the diagonals of a matrix
    deg = np.asarray(P.sum(axis=0)).flatten()
    deg[np.isinf(deg)] = 0
    D = diags(deg)
    L = D - P
    return L


def runGeneMANIA(L, y):
    """
    *L*: Laplacian of the original network
    *y*: vector of positive and negative assignments
    """
    y = y.copy()
    # setup the y vector with the value for the unknowns
    num_pos = len(np.where(y == 1)[0])
    num_neg = len(np.where(y == -1)[0])
    # taken from the GeneMANIA paper
    k = (num_pos - num_neg) / float(num_pos + num_neg)
    y[np.where(y == 0)[0]] = k

    start_time = time.process_time()
    # try solving for s directly
    M = eye(L.shape[0]) + L
    # use scipy's conjugate gradient solver
    f, info = cg(M, y)
    #f = spsolve(M, y)
    #info = ''
    total_time = time.process_time() - start_time
    #print("Solved GeneMANIA using conjugate gradient (%0.2f sec). Info: %s, k=%0.2f" % (
    #    total_time, str(info), k))

    return f, total_time
