# Python implementation of GeneMANIA

import time
#from tqdm import trange, tqdm
import alg_utils
import numpy as np
from scipy.sparse import eye, diags
from scipy.sparse.linalg import spsolve, cg, spilu, LinearOperator
from scipy.io import savemat
import os
import matlab.engine


def setup_laplacian(W):
    # first normalize the network
    P = alg_utils._net_normalize(W)
    # take the column sum and set them as the diagonals of a matrix
    deg = np.asarray(P.sum(axis=0)).flatten()
    deg[np.isinf(deg)] = 0
    D = diags(deg)
    L = D - P
    return L


def runGeneMANIA(eng, L, y, tol=1e-05, verbose=False):
    """
    *L*: Laplacian of the original network
    *y*: vector of positive and negative assignments
    *tol*: Conjugate Gradient tolerance for convergance
    """
    y_orig = y
    y = y.copy()
    # setup the y vector with the value for the unknowns
    num_pos = len(np.where(y == 1)[0])
    num_neg = len(np.where(y == -1)[0])
    # taken from the GeneMANIA paper
    k = (num_pos - num_neg) / float(num_pos + num_neg)
    y[np.where(y == 0)[0]] = k
    #M = eye(L.shape[0]) + L
    eng.eval("M = speye(size(L)) + L;", nargout=0)

    eng.workspace['y'] = matlab.double(list(y))
    eng.workspace['y_orig'] = matlab.double(list(y_orig))
    eng.workspace['tol'] = tol
    # need to convert to transpose to get the right dimensions
    eng.eval("y = y';", nargout=0)
    eng.eval("y_orig = y_orig';", nargout=0)


    #start_process_time = time.process_time()
    start_wall_time = time.time()
    eng.eval("tic; t=cputime; [f,FLAG,RELRES,ITER] = pcg(M, y); process_time = cputime - t; wall_time = toc;", nargout=0)
    #eng.eval("tic; t=cputime; [f,FLAG,RELRES,ITER] = cgs(M, y); process_time = cputime - t; wall_time = toc;", nargout=0)
    # also try the original GM matlab code
    #eng.eval("addpath('src/algorithms/weight_networks/GeneMANIACode_v2/');")
    #eng.eval("tic; t=cputime; [f,RELRES,ITER] = getScoreVectorCG(y_orig, M); process_time = cputime - t; wall_time = toc;", nargout=0)
    #eng.eval("tic; t=cputime; finit=zeros(size(y,1),1); maxit=1000; [f,RELRES,ITER] = conjGrad(finit,y,M,maxit,tol); RELRES=sum(RELRES.^2); process_time = cputime - t; wall_time = toc;", nargout=0)

    #flag = eng.workspace['FLAG']
    relres = eng.workspace['RELRES']
    iters = eng.workspace['ITER']
    wall_time = eng.workspace['wall_time']
    process_time = eng.workspace['process_time']

    #f, info = cg(M, y, tol=tol)
    #process_time = time.process_time() - start_process_time
    py_wall_time = time.time() - start_wall_time
    if verbose:
        print("Solved GeneMANIA using pcg. RELRES: %s, ITER: %s, wall_time: %s, process_time: %s, py_wall_time: %s" % (
            relres, iters, wall_time, process_time, py_wall_time))
        #print("Solved GeneMANIA using pcg. FLAG: %s, RELRES: %s, ITER: %s, wall_time: %s, py_wall_time: %s, py_process_time: %s" % (
        #    flag, relres, iters, wall_time, py_wall_time, process_time))
    f = np.asarray(eng.workspace["f"]._data)

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

    return f, process_time, wall_time, iters
