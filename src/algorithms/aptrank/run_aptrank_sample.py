import numpy as np
from scipy.io import loadmat, savemat
import scipy.sparse as sparse
import scipy as sp
from aptrank import AptRank


# The equivalent function for splitRT.m
def shuffleSplit(r, t):

    if t < 0 or t > 1:
        print("T must be in [0,1]")
    m, n = r.shape
    ri, rj, rv = sparse.find(r)
    rlen = len(ri)
    p = np.random.permutation(rlen)
    p = np.delete(p, np.arange(0, rlen - int(t*rlen)))
    cp = np.setdiff1d(np.arange(0, rlen), p)

    # Rtrain = sparse.coo_matrix((ri[p], rj[p], rv[p]), shape=(m, n), dtype=float).tocsr()
    Rtrain = sparse.csr_matrix((rv[p], (ri[p], rj[p])), shape=(m, n), dtype=float)
    Rtest = sparse.csr_matrix((rv[cp], (ri[cp], rj[cp])), shape=(m, n), dtype=float)


    # spliter = ShuffleSplit(n_splits=1, test_size=self.T, random_state=None)
    # indices = list(spliter.split(self.Rtrain))
    # train_indices = indices[0]
    # test_indices = indices[1]

    return Rtrain, Rtest

# lamb percent of H going to the bottom
# the rest percent of H^T point to the top
def directH(H, lamb):

    m, n = H.shape
    if m != n:
        print("H must be squared")
        exit(2)
    dH = colstonnz(lamb*H) + colstonnz((1 - lamb)*H.T)
    return dH


"""
Convert matrix B into a column stochastic matrix
Normalize each column by the number of non-zero entries in that column
"""
def colstonnz(B):

    C = B.copy()
    if not sparse.isspmatrix(C):
        C = sparse.csr_matrix(C)
    m,n = C.shape
    bi, bj, bv = sparse.find(C)
    d = sp.zeros(n)
    for i in np.arange(0, n):
        # counting amount of nonzero entries in each row
        d[i] = len(sparse.find(B[:, i])[0])

    divsor = d[bj]
    val = np.divide(bv, divsor)

    C = sparse.csr_matrix((val, (bi, bj)), shape=(m, n), dtype=float)
    return C


def newROC(targets, outputs):
    # print(outputs.shape)
    numClasses, numSamples = targets.shape
    numClasses2, numSamples2 = outputs.shape

    if numClasses != numClasses2 or numSamples != numSamples2:
        print("Targets and outputs have different dimenstions")
        exit(4)

    if numClasses == 1:
        targets = np.vstack((targets.toarray(), np.ones(numSamples) - targets.toarray()))
        outputs = np.vstack((outputs, 1 - outputs - np.finfo(float).eps * (outputs == 0.5)))
        tpr, fpr, thresholds = newROC(targets, outputs)
        tpr = tpr[0]
        fpr = fpr[0]
        thresholds = thresholds[0]

        return tpr, fpr, thresholds

    fpr = np.empty(numClasses, dtype=object)
    tpr = np.empty(numClasses, dtype=object)
    thresholds = np.empty(numClasses, dtype=object)

    for i in np.arange(0, numClasses):
        print("starting the %dth class" % i)
        tpr[i], fpr[i], thresholds[i] = roc_one(targets[i, :], outputs[i, :])
        print("completed the %dth class" % i)

    return tpr, fpr, thresholds


"""
Helper function for NewROC
"""
def roc_one(targets, outputs):

    numSamples = len(targets)
    numPosTargets = np.sum(targets)
    numNegTargets = numSamples - numPosTargets

    thresholds = np.unique(np.concatenate([[0], outputs, [1]]))
    numThresholds = len(thresholds)

    sortedPosTargetoutputs = np.sort(outputs[targets == 1])
    numPosTargetOutputs = len(sortedPosTargetoutputs)
    sortedNegTargetOutputs = np.sort(outputs[targets == 0])
    numNegTargetsOutputs = len(sortedNegTargetOutputs)

    fpcount = np.zeros((1, numThresholds))
    tpcount = np.zeros((1, numThresholds))

    posInd = 0
    negInd = 0
    print("roc_one")
    for i in np.arange(numThresholds):
        threshold = thresholds[i]
        while posInd < numPosTargetOutputs and sortedPosTargetoutputs[posInd] <= threshold:
            posInd += 1
        tpcount[0, i] = numPosTargetOutputs - posInd + 1
        while negInd < numNegTargetsOutputs and sortedNegTargetOutputs[negInd] <= threshold:
            negInd += 1
        fpcount[0, i] = numNegTargetsOutputs - negInd + 1
    tpr = np.fliplr(tpcount) / numPosTargets
    fpr = np.fliplr(fpcount) / numNegTargets
    print("roc_one done")
    return tpr, fpr, thresholds


"""
Python version calAUC
Test set evaluation without using the positive samples in Rtrain.
No figure display, AUROC calculated only.
The input 3 matrices must have the same dimensions.
"""


def calcAUC(X, Rtrain, Rtest):
    # remove empty rows and columns in Rtest
    cx = sparse.find(Rtest.sum(axis=0))[1]
    rx = sparse.find(Rtest.sum(axis=1))[0]
    Rtest = Rtest[np.ix_(rx, cx)]
    Rtrain = Rtrain[np.ix_(rx, cx)]
    X = X[np.ix_(rx, cx)]
    Rtrain = np.reshape(Rtrain, (Rtrain.shape[0] * Rtrain.shape[1], 1), order='F')
    Rtest = np.reshape(Rtest, (Rtest.shape[0] * Rtest.shape[1], 1), order='F')
    X = np.reshape(X, (X.shape[0] * X.shape[1], 1), order='F')
    del_ind = sparse.find(Rtrain)[0]
    # del_ind = np.setdiff1d(ind_nnz, np.arange(Rtrain.shape[0]))
    Rtest.data[np.where(np.in1d(Rtest.col, del_ind))] = 0
    Rtest.eliminate_zeros()
    X[del_ind] = 0
    X = X + np.random.rand(X.shape[0], X.shape[1]) * np.finfo(float).eps  # avoid equal rankings
    tpr, fpr, c = newROC(Rtest.T, X.T)
    auc = 0

    for i in np.arange(len(c) - 1):
        auc = auc + (fpr[0, i + 1] - fpr[0, i]) * (tpr[0, i] + tpr[0, i + 1]) / 2
    return auc

"""
Python version calMAP
CAlculate Mean Average Precision
Do not count the postive entries in Rtrain
The input 3 matrices must have the same dimensions.
"""
def calcMAP(X, Rtrain, Rtest):

    cx = sparse.find(Rtest.sum(axis=0))[1]
    rx = sparse.find(Rtest.sum(axis=1))[0]
    Rtest = Rtest[np.ix_(rx, cx)]
    X = X[np.ix_(rx, cx)]
    if Rtrain.size > 0:
        Rtrain = Rtrain[np.ix_(rx, cx)]
    m = X.shape[0]
    ap = np.zeros((m, 1))
    Q = []
    for i in np.arange(0, m):
        if Rtrain.size > 0:
            rn0 = sparse.find(Rtrain[i, :] == 0)[1]
            Q = np.vstack((X[i, rn0], Rtest[i, rn0].toarray()))
        else:
            Q = np.vstack(([X[i, :], Rtest[i, :].toarray()]))
        Q = Q.T
        Q = Q[Q[:, 0].argsort(), ][::-1]
        tp = sparse.find(Q[:, 1])[1]
        if len(tp) > 0:
            u = np.arange(1, len(tp) + 1)
            u = u.T
            prec = u/(tp + 1)
            ap[i] = np.mean(prec)
        else:
            ap[i] = 0
    return np.mean(ap, dtype=float)


rho = 0.5
data = loadmat("2018_06-seq-sim-e1e-25-aptrank-bp-annotations-and-go-dag.mat")
# data = np.load("2018_06-seq-sim-e1e-25-net.npz")
# data_sets = loadmat("checker2.mat")
# data_sets = loadmat("init.mat")
G = data['G']
H = data['H']
R = data['R']
# Rtrain = data_sets['Rtrain']
# Rtest = data_sets['Rtest']
Rtrain, Rtest = shuffleSplit(R, rho)

lam = 0.5
dh = directH(H, lam)
K = 2
S = 1
T = 0.5
diffusion_type = 'twoway'
apt = AptRank(G, Rtrain.T, dh, K, S, T, 12,diffusion_type)
xa = apt.algorithm()
# results = loadmat("result_test.mat")['Xa']
# diff = xa - results
# print("Max difference is %0.10f" % diff.max())
# print("Difference mean is %0.10f" % diff.mean())
# print("Difference std is %0.10f" % diff.std())
# savemat("difference.mat", diff)
print("AUC = %0.7f" % calcAUC(xa, Rtrain.T, Rtest.T))
print("Map = %0.7f" % calcMAP(xa, Rtrain.T, Rtest.T))
