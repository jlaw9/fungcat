# Python implementation of SinkSource

import time
#from tqdm import trange, tqdm
import alg_utils
import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.sparse import linalg
from scipy.io import savemat
import os
import matlab.engine


def SinkSource_matlab(eng, P, f, max_iters=1000, eps=0.0001, a=0.8, verbose=False):
    """ Iterate over the graph until all of the scores of the unknowns have converged
    *max_iters*: maximum number of iterations
    *eps*: Epsilon convergence cutoff. If the maximum node score change from one iteration to the next is < *eps*, then stop
    *a*: alpha is the insulation parameter

    *returns*: scores array, process_time, wall_time, 
        # of iterations needed to converge, total # of computations performed
    """
    # total number of computations performed during the update function
    total_comp = 0
    # f technically already contains the first iteration
    #s = f.copy()
    # UPDATE 2018-12-21: just start at zero so the time taken measured here is correct
    s = np.zeros(len(f))
    prev_s = s.copy()
    # convert these to matlab variables
    eng.workspace['f'] = matlab.double(list(f))
    eng.workspace['s'] = matlab.double(list(s))
    eng.workspace['prev_s'] = matlab.double(list(prev_s))
    # need to convert to transpose to get the right dimensions
    eng.eval("f = f'; s = s'; prev_s = prev_s';", nargout=0)
    eng.workspace['a'] = a
    eng.workspace['ss_eps'] = eps
    eng.workspace['max_iters'] = max_iters

    # now run SinkSource
    # solve using power iteration
    # for some reason, taking the transpose here and then again before calling the dot product is about 3x faster
    eng.eval("M = a*M';", nargout=0)
    ss_matlab = """tic; for iters = 1:max_iters
    s = M'*prev_s + f;
    max_d = max(s - prev_s);
    if max_d <= ss_eps
       break
    end
    prev_s = s;
end
wall_time = toc;
"""
    wall_time = time.time()
    eng.eval(ss_matlab, nargout=0)

#   # solve the equation directly
#   # rewrite as Ax = b
#   # P' is included to undo the transpose above
#   eng.eval("I = speye(size(P));", nargout=0)
#   eng.eval("A = (I - a*P);", nargout=0)
#   #ss_matlab = "tic; s = A\f; wall_time = toc; iters = -1;"
#   wall_time = time.time()
#   ss_matlab = "tic; s = mldivide(A,f); wall_time = toc; iters = -1;"
#   eng.eval(ss_matlab, nargout=0)

    #iters = eng.eval(ss_matlab)
#    for iters in range(1,max_iters+1):
#        update_time = time.time()
#        #s = a*csr_matrix.dot(P, prev_s) + f
#        eng.eval("s = a*P*prev_s + f;", nargout=0)
#
#        # TODO store this in a list and return it
#        #max_d = (s - prev_s).max()
#        max_d = eng.eval("max(s - prev_s);")
#        if verbose:
#            print("\t\titer %d; %0.4f sec to update scores, max score change: %0.6e" % (iters, time.time() - update_time, max_d))
#        if max_d <= eps:
#            # converged!
#            break
#        #prev_s = s.copy()
#        eng.eval("prev_s = s;", nargout=0)

    wall_time2 = time.time() - wall_time
    wall_time = eng.workspace['wall_time']
    iters = eng.workspace['iters']
    #process_time = time.process_time() - process_time
    if verbose:
        #print("SinkSource converged after %d iterations (%0.3f sec, %0.3f process time)" % (iters, wall_time, process_time))
        print("SinkSource converged after %d iterations (%0.3f sec, %0.3f python sec)" % (iters, wall_time, wall_time2))
    s = np.asarray(eng.workspace["s"]._data)

    return s, wall_time, wall_time, iters


def runSinkSource(eng, P, positives, negatives=None, max_iters=1000, eps=0.0001, a=0.8, verbose=False):
    """
    *P*: Network as a scipy sparse matrix. Should already be normalized
    *positives*: numpy array of node ids to be used as positives
    *negeatives*: numpy array of node ids to be used as negatives. 
        If not given, will be run as SinkSourcePlus. 
        For SinkSourcePlus, Lambda should already be included in the graph normalization process. 
        See the function normalizeGraphEdgeWeights in alg_utils.py 
    *max_iters*: max # of iterations to run SinkSource. 
        If 0, use spsolve to solve the equation directly 
    """
    num_nodes = P.shape[0]
    # TODO check to make sure the graph is normalized 
    # remove the positive and negative nodes from the graph 
    # and setup the f vector which contains the influence from positive and negative nodes
    f, node2idx, idx2node = setupScores(
        eng, P, positives, negatives, a=a, remove_nonreachable=False)
    comp = 0
    ## make sure we're getting enough precision
    #f = f.astype('float128')

    #if max_iters > 0:
    s, process_time, wall_time, num_iters = SinkSource_matlab(
        eng, P, f, max_iters=max_iters, eps=eps, a=a, verbose=verbose)
#    else:
#        # Solve for s directly. Scipy uses the form Ax=b to solve for x
#        # SinkSource equation: (I - aP)s = f
#        # TODO add an option to specify which solver to use
#        # spsolve basically stalls for denser networks (e-value cutoff 0.1)
#        #solver = "spsolve"
#        solver = "cg"
#        # eye is the identity matrix
#        A = eye(P.shape[0]) - a*P
#        #print("testing if the matrix is positive definite. 10 largest and 10 smallest eigenvalues:")
#        #print(linalg.eigs(A, k=10, which="LM", return_eigenvectors=False))
#        #print(linalg.eigs(A, k=10, which="SM", return_eigenvectors=False))
#        # this measures the amount of time taken by all processors
#        # and Conjugate Gradient is paralelized
#        start_process_time = time.process_time()
#        # this measures the amount of time that has passed
#        start_wall_time = time.time()
#        if solver == 'spsolve':
#            num_iters = 0
#            s = linalg.spsolve(A, f)
#        elif solver == 'cg':
#            # keep track of the number of iterations
#            num_iters = 0
#            def callback(xk):
#                #nonlocal num_iters
#                num_iters += 1
#            # cg (conjugate gradient) is the fastest of the other available solvers
#            s, info = linalg.cg(A, f, callback=callback, maxiter=1000)
#        process_time = time.process_time() - start_process_time 
#        wall_time = time.time() - start_wall_time
#        if verbose:
#            print("Solved SS using %s (%0.3f sec, %0.3f process time). %s" % (
#                solver, wall_time, process_time, 
#                "iters: %d, max_iters: %d" % (num_iters, 1000) if solver == 'cg' else ''))
#        comp = 0

    # map back from the indices after deleting pos/neg to the original indices
    scores_arr = np.zeros(num_nodes)
    indices = [idx2node[n] for n in range(len(s))]
    scores_arr[indices] = s

    return scores_arr, process_time, wall_time, num_iters, comp


def runLocal(P, positives, negatives=None):
    f = np.zeros(P.shape[0])
    f[positives] = 1
    if negatives is not None:
        f[negatives] = -1

    start_process_time = time.process_time()
    s = csr_matrix.dot(P, f)
    process_time = time.process_time() - start_process_time

    return s, process_time


def setupScores(eng, P, positives, negatives=None, a=1, remove_nonreachable=True, verbose=False):
    """
    Updated to remove nodes from the P matrix inside of the matlab engine
    """
    #print("Initializing scores and setting up network")
    pos_vec = np.zeros(P.shape[0])
    pos_vec[positives] = 1
    #if negatives is not None:
    #    pos_vec[negatives] = -1

    # f contains the fixed amount of score coming from positive nodes
    f = a*csr_matrix.dot(P, pos_vec)

    fixed_nodes = positives 
    if negatives is not None:
        fixed_nodes = np.concatenate([positives, negatives])
    node2idx, idx2node = alg_utils.build_index_map(range(len(f)), set(list(fixed_nodes)))
    # removing the fixed nodes is slightly faster than selecting the unknown rows
    # remove the fixed nodes from the graph
    fixed_nodes = np.asarray(list(fixed_nodes)) if not isinstance(fixed_nodes, np.ndarray) else fixed_nodes
    eng.workspace['fixed_nodes'] = matlab.double(list(fixed_nodes))
    # convert to matlab indices
    eng.eval("fixed_nodes = fixed_nodes'; fixed_nodes = fixed_nodes + 1;", nargout=0)
    # now remove the nodes from the P matrix
    eng.eval("M = P; M(fixed_nodes,:) = []; M(:,fixed_nodes) = [];", nargout=0)
    #P = delete_nodes(P, fixed_nodes)
    # and from f
    f = np.delete(f, fixed_nodes)
    #assert P.shape[0] == P.shape[1], "Matrix is not square"
    #assert P.shape[1] == len(f), "f doesn't match size of P"

    return f, node2idx, idx2node


