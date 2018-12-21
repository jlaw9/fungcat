# Python implementation of SinkSource

import time
#from tqdm import trange, tqdm
import alg_utils
import numpy as np
from scipy.sparse import csr_matrix, eye
from scipy.sparse import linalg
#try:
#    # this only worked installed from conda 
#    from numba import njit, prange
#    use_numba = True
#    print("\tUsing numba for parallel matrix multiplication")
#except ImportError:
#    use_numba = False
# don't use numba for now
use_numba = False


if use_numba:
    # This was about 2.5x faster than regular sinksource dot product in my tests
    # due to the parallelization.
    # Results are mixed in practice. With a large 70M edges network, I got ~2x speedup
    # copied from here: https://stackoverflow.com/questions/46924092/how-to-parallelize-this-python-for-loop-when-using-numba
    @njit(parallel=True)
    def csr_dot_numba(x, Adata, Aindices, Aindptr, Ashape):

        numRowsA = Ashape[0]    
        Ax       = np.zeros(numRowsA)

        for i in prange(numRowsA):
            Ax_i = 0.0    
            for dataIdx in np.arange(Aindptr[i],Aindptr[i+1]):

                j     = Aindices[dataIdx]
                Ax_i += Adata[dataIdx] * x[j]

            Ax[i] = Ax_i                

        return Ax


# old option: *scores*: rather than start from scratch, we can start from a previous run's (or different GO term's) scores
def SinkSource_scipy(P, f, max_iters=1000, eps=0.0001, a=0.8, verbose=False):
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

    wall_time = time.time()
    # this measures the amount of time taken by all processors
    process_time = time.process_time()
    for iters in range(1,max_iters+1):
        update_time = time.time()
        if use_numba:
            s = a*csr_dot_numba(prev_s, P.data, P.indices, P.indptr, P.shape) + f
        else:
            s = a*csr_matrix.dot(P, prev_s) + f

        # TODO store this in a list and return it
        max_d = (s - prev_s).max()
        if verbose:
            print("\t\titer %d; %0.4f sec to update scores, max score change: %0.6e" % (iters, time.time() - update_time, max_d))
        if max_d <= eps:
            # converged!
            break
        prev_s = s.copy()

    wall_time = time.time() - wall_time
    process_time = time.process_time() - process_time
    total_comp += len(P.data)*iters
    if verbose:
        print("SinkSource converged after %d iterations (%0.3f sec, %0.3f process time)" % (iters, wall_time, process_time))

    return s, process_time, wall_time, iters, total_comp


def runSinkSource(P, positives, negatives=None, max_iters=1000, eps=0.0001, a=0.8, verbose=False):
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
    P, f, node2idx, idx2node = alg_utils.setupScores(
        P, positives, negatives, a=a, remove_nonreachable=False)
    ## make sure we're getting enough precision
    #f = f.astype('float128')

    if max_iters > 0:
        s, process_time, wall_time, num_iters, comp = SinkSource_scipy(
            P, f, max_iters=max_iters, eps=eps, a=a, verbose=verbose)
    else:
        # Solve for s directly. Scipy uses the form Ax=b to solve for x
        # SinkSource equation: (I - aP)s = f
        # TODO add an option to specify which solver to use
        # spsolve basically stalls for denser networks (e-value cutoff 0.1)
        #solver = "spsolve"
        solver = "cg"
        # eye is the identity matrix
        A = eye(P.shape[0]) - a*P
        #print("testing if the matrix is positive definite. 10 largest and 10 smallest eigenvalues:")
        #print(linalg.eigs(A, k=10, which="LM", return_eigenvectors=False))
        #print(linalg.eigs(A, k=10, which="SM", return_eigenvectors=False))
        # this measures the amount of time taken by all processors
        # and Conjugate Gradient is paralelized
        start_process_time = time.process_time()
        # this measures the amount of time that has passed
        start_wall_time = time.time()
        if solver == 'spsolve':
            num_iters = 0
            s = linalg.spsolve(A, f)
        elif solver == 'cg':
            # keep track of the number of iterations
            num_iters = 0
            def callback(xk):
                nonlocal num_iters
                num_iters += 1
            # cg (conjugate gradient) is the fastest of the other available solvers
            s, info = linalg.cg(A, f, callback=callback, maxiter=1000)
        process_time = time.process_time() - start_process_time 
        wall_time = time.time() - start_wall_time
        if verbose:
            print("Solved SS using %s (%0.3f sec, %0.3f process time). %s" % (
                solver, wall_time, process_time, 
                "iters: %d, max_iters: %d" % (num_iters, 1000) if solver == 'cg' else ''))
        comp = 0

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
