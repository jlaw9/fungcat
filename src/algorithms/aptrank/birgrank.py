
#import numpy as np
import algorithms.gain_scipy.alg_utils as alg_utils
from scipy.sparse import csc_matrix, vstack, hstack, eye
from scipy.sparse.linalg import spsolve
# try using a conjugate gradient solver 
# unfortunately, gives an "incompatible dimensions" error
# likely because it's not a dense matrix
from scipy.sparse.linalg import cg
import time
#import sys
#from tqdm import tqdm, trange

__author__ = "Jeff Law"

# python implementation of the BirgRank algorithm 
# see https://github.rcac.purdue.edu/mgribsko/aptrank for the original matlab implementation
# full citation:
#Jiang, B., Kloster, K., Gleich, D. F., & Gribskov, M., Aptrank: an adaptive pagerank model for protein function prediction on bi-relational graphs, Bioinformatics, 33(12), 1829-1836 (2017).  http://dx.doi.org/10.1093/bioinformatics/btx029


def birgRank(G, Rtrain, dH, alpha=.5, theta=.5, mu=.5):
    """
    *Rtrain*: Matrix containing 0s and 1s for annotations. 
        Rows are nodes, columns are goids
    *alpha*: teleportation parameter in Personalized PageRank.
    *theta*: (1-theta) percent of Rtrain used in seeding vectors.
    *mu*: (1-mu) percent of random walkers diffuse from G via Rtrain to H.

    *returns*: Xh - m*n matrix of scores for each GO term, protein pair, same size as Rtrain
    """

    #G = G.tocsc()
    #Rtrain = Rtrain.tocsc()
    #dH = dH.tocsc()
    print("Starting birgRank")
    m, n = Rtrain.shape

    # make seeding vectors (matrix)
    B = vstack([theta*eye(m), (1-theta)*Rtrain.T]).tocsc()
    # column normalize the matrix
    B = alg_utils.normalizeGraphEdgeWeights(B, axis=0)
    B = (1-alpha)*B

    # make transition matrix
    G = G/2 + G.T/2  # make sure G is symmetric
    # final matrix is:
    # [G  0]
    # [RT H]
    P = vstack([hstack([mu*G, csc_matrix((m,n))]), hstack([(1-mu)*Rtrain.T, dH])]).tocsc()
    # normalize using the same normalization as AptRank
    P = alg_utils.normalizeGraphEdgeWeights(P, axis=0)  # column normalization
    A = eye(m+n) - alpha*P

    print("Starting solver")
    start_time = time.process_time()

    print("solving for Xg")
    # make sure they're in csc format for the solvers
    A = A.tocsc()
    B = B.tocsc()
    # TODO solve PageRank linear system using 3-block solver
    # now solve the linear system X = A/B
    # split it up into two equations to solve
    # (I-alpha*G)Xg = (1-alpha)I
    Xg = spsolve(A[:m,:][:,:m], B[:m,:])
    # gives a memory error
    #Xg = cg(A[:m,:][:,:m], B[:m,:].todense())
    print("solving for Xh")
    # alpha*RT*Xg = (I - alpha*H)Xh
    Xh = spsolve(A[m:,:][:,m:], (B[m:,:] - A[m:,:][:,:m]*Xg))
    #Xh = cg(A[m:,:][:,m:], (B[m:,:] - A[m:,:][:,:m]*Xg).todense())

    total_time = time.process_time() - start_time
    print("Solved birgRank using sparse linear system (%0.2f sec)" % (total_time))

    # transpose Xh so it has the same dimensions as Rtrain
    return Xh.T
