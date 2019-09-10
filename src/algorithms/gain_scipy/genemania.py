# Python implementation of GeneMANIA

import time
#from tqdm import trange, tqdm
import alg_utils
import numpy as np
from scipy.sparse import eye, diags
from scipy.sparse.linalg import spsolve, cg, spilu, LinearOperator


def setup_laplacian(W):
    # first normalize the network
    P = alg_utils._net_normalize(W)
    # take the column sum and set them as the diagonals of a matrix
    deg = np.asarray(P.sum(axis=0)).flatten()
    deg[np.isinf(deg)] = 0
    D = diags(deg)
    L = D - P
    return L


def runGeneMANIA(L, y, tol=1e-05, verbose=False):
    """
    *L*: Laplacian of the original network
    *y*: vector of positive and negative assignments
    *tol*: Conjugate Gradient tolerance for convergance
    """
    y = y.copy()
    # setup the y vector with the value for the unknowns
    num_pos = len(np.where(y == 1)[0])
    num_neg = len(np.where(y == -1)[0])
    # taken from the GeneMANIA paper
    k = (num_pos - num_neg) / float(num_pos + num_neg)
    y[np.where(y == 0)[0]] = k

    start_process_time = time.process_time()
    start_wall_time = time.time()
    M = eye(L.shape[0]) + L
    # use scipy's conjugate gradient solver
    # keep track of the number of iterations
    num_iters = 0
    # also keep track of the maximum difference of scores from one iteration to the next (epsilon) like we do for SinkSource
    #max_d = 0
    #prev_xk = np.zeros(len(y))
    def callback(xk):
        nonlocal num_iters
        num_iters += 1
        #nonlocal max_d
        #nonlocal prev_xk
        #max_d = abs(xk - prev_xk).max()
        #prev_xk = xk.copy() 

    f, info = cg(M, y, tol=tol, callback=callback)
    #f = spsolve(M, y)
    #info = ''
    process_time = time.process_time() - start_process_time
    wall_time = time.time() - start_wall_time
    if verbose:
        #print("Solved GeneMANIA using conjugate gradient (%0.2f sec, %0.2f sec process_time). Info: %s, k=%0.2f, max_d=%0.4e, iters=%d" % (
        #    wall_time, total_time, str(info), k, max_d, num_iters))
        print("Solved GeneMANIA using conjugate gradient (%0.2f sec, %0.2f sec process_time). Info: %s, k=%0.2f, iters=%d" % (
            wall_time, process_time, str(info), k, num_iters))

    # TODO try to pre-process the M matrix to be able to solve for multiple GO terms faster
    # using the pre-process matrix isn't actually faster according to the tests I've done so far
#
#    M = M.tocsc()
#    start_time = time.process_time()
#    s_time = time.time()
#    M2 = spilu(M)
#    M2 = LinearOperator(M.shape, M2.solve)
#    print("Incomplete LU decomposition: %0.2f sec (%0.2f process_time)" % (time.time() - s_time, time.process_time() - start_time))
#
#    start_time = time.process_time()
#    f, info = cg(M, y, M=M2)
#    print("Solved GeneMANIA with approximate inverse of A: %0.2f sec (%0.2f process_time)" % (time.time() - s_time, time.process_time() - start_time))

    return f, process_time, wall_time, num_iters
