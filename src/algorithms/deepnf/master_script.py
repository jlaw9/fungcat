from optparse import OptionParser
import os
import numpy as np
os.environ["KERAS_BACKEND"] = "tensorflow"
from keras.models import load_model, Model
import scipy.sparse as sp
from scipy.io import savemat, loadmat
import sys
from sklearn.preprocessing import minmax_scale
from deepNF_class import DeepNF
import matplotlib.pyplot as plt


def parse_args(args):

    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    # general parameters
    parser.add_option("", "--models_path",
                      help="model path")
    parser.add_option("", "--results_path",
                      help="results path")
    parser.add_option("", "--annot",
                      help="annotation file")
    parser.add_option("", "--input_networks",
                      help="input network files in the form of an array, seperated by comma."
                           "--input_networks=/data/file1.mat,/data/file2.mat")
    parser.add_option("", "--prefix",
                      help="prefix")
    parser.add_option("", "--goids",
                      help="file containing go ids")
    parser.add_option("", "--protids",
                      help="file containing protein ids")

    # parameters with default value
    parser.add_option("", "--valid_type", action="store_true",
                      help="Validation method: Temporal holdout or Crossvalidation", default=True)
    parser.add_option("", "--epochs", type='int',
                      help="iterations for trainning deep learning model", default=10)
    parser.add_option("", "--batch_size", type='int',
                      help="sample size for trainning each iteration", default=128)
    parser.add_option("", "--ntrials", type='int',
                      help="number of cv trials", default=10)
    parser.add_option("", "--alpha", type="int",
                      help="propagation parameter", default=0.98)
    parser.add_option("", "--select_arch",
                      help="a number 1-6 (see below)")

    # optional parameters
    parser.add_option("-k", "--topK", type='int',
                      help="top K prediction score for protein go term pair")
    parser.add_option("", "--cutoff_hi", type='int',
                      help="upper bound for cutoff")
    parser.add_option("", "--cutoff_lo", type='int',
                      help="lower bound for cutoff")
    parser.add_option("", "--prep", action="store_true", default=False,
                      help="Option to preprocess the input network files")

    (options, args) = parser.parse_args(args)

    if options.valid_type is None or not "cv" or not "th":
        print "invalid validation type: %s" % options.valid_type
        sys.exit(1)

    if options.models_path is None or not os.path.isdir(options.models_path):
        print "please correct path for models, either path is not set or output directory doesn't exist"
        sys.exit(2)

    if options.results_path is None or not os.path.isdir(options.results_path):
        print "please correct path for results, either path is not set or output directory doesn't exist"
        sys.exit(3)

    if options.annot is None or not os.path.exists(options.annot):
        print "please provide the annotation file/" \
              "the correct path to the annotation file in order for deepNF to validate"
        sys.exit(4)

    if options.goids is None or not os.path.exists(options.goids):
        print "please correct path for go ids, either path is not set or output directory doesn't exist"
        sys.exit(5)

    if options.protids is None or not os.path.exists(options.protids):
        print "please correct path for protein ids, either path is not set or output directory doesn't exist"
        sys.exit(6)

    if options.input_networks is None:
        print "please provide the network files"
        sys.exit(7)

    nets = options.input_networks.split(",")
    for net in nets:
        if not os.path.exists(net):
            print "network file %s doesn't exist, please check the file" % net
            sys.exit(8)

    return options


def _load_network(input, mtrx='adj'):
    print("### Loading [%s]..." % ("current"))
    if mtrx == 'adj':
        # i, j, val = np.loadtxt(filename).T
        # A = coo_matrix((val, (i-1, j-1)))
        # id are starts at 0, so no need to decrease the number
        # A = coo_matrix((val.astype(int), (i.astype(int), j.astype(int))))
        # A = A.todense()
        # A = np.squeeze(np.asarray(A))
        A = input.toarray()
        if A.min() < 0:
            print("### Negative entries in the matrix are not allowed!")
            A[A < 0] = 0
            print("### Matrix converted to nonnegative matrix.")
            print
        if np.allclose(A.T, A, atol=1e-8):
            # if A.T == A:
            pass
        else:
            print("### Matrix not symmetric!")
            A = A + A.T
            print("### Matrix converted to symmetric.")
    else:
        print("### Wrong mtrx type. Possible: {'adj', 'inc'}")
    A = A - np.diag(np.diag(A))
    A = A + np.diag(A.sum(axis=1) == 0)

    return A


def load_networks(input_net, mtrx='adj'):
    """
    Function for loading Mashup files
    Files can be downloaded from:
        http://cb.csail.mit.edu/cb/mashup/
    """
    Nets = []
    for net in input_net:
        Nets.append(_load_network(net, mtrx))

    return Nets


def _net_normalize(X):
    """
    Normalizing networks according to node degrees.
    """
    if X.min() < 0:
        print("### Negative entries in the matrix are not allowed!")
        X[X < 0] = 0
        print("### Matrix converted to nonnegative matrix.")
        print()
    if (X.T == X).all():
        pass
    else:
        print("### Matrix not symmetric.")
        X = X + X.T - np.diag(np.diag(X))
        print("### Matrix converted to symmetric.")

    # normalizing the matrix
    deg = X.sum(axis=1).flatten()
    deg = np.divide(1., np.sqrt(deg))
    deg[np.isinf(deg)] = 0
    D = np.diag(deg)
    X = D.dot(X.dot(D))

    return X


def net_normalize(Net):
    """
    Normalize Nets or list of Nets.
    """
    if isinstance(Net, list):
        for i in range(len(Net)):
            Net[i] = _net_normalize(Net[i])
    else:
        Net = _net_normalize(Net)

    return Net


def _scaleSimMat(A):
    """Scale rows of similarity matrix"""
    A = A - np.diag(np.diag(A))
    A = A + np.diag(A.sum(axis=0) == 0)
    col = A.sum(axis=0)
    A = A.astype(np.float)/col[:, None]

    return A


def RWR(A, K=3, alpha=0.98):
    """Random Walk on graph"""
    A = _scaleSimMat(A)
    # Random surfing
    n = A.shape[0]
    P0 = np.eye(n, dtype=float)
    P = P0.copy()
    M = np.zeros((n, n), dtype=float)
    for i in range(0, K):
        P = alpha*np.dot(P, A) + (1. - alpha)*P0
        M = M + P

    return M


def PPMI_matrix(M):
    """ Compute Positive Pointwise Mutual Information Matrix"""
    M = _scaleSimMat(M)
    n = M.shape[0]
    col = np.asarray(M.sum(axis=0), dtype=float)
    col = col.reshape((1, n))
    row = np.asarray(M.sum(axis=1), dtype=float)
    row = row.reshape((n, 1))
    D = np.sum(col)

    np.seterr(all='ignore')
    PPMI = np.log(np.divide(D*M, np.dot(row, col)))
    PPMI[np.isnan(PPMI)] = 0
    PPMI[PPMI < 0] = 0

    return PPMI


def run():

    opt = parse_args(sys.argv)

    # reading annotation file
    annot = sp.load_npz(opt.annot).transpose().toarray()
    # TODO teneray classification
    annot[annot < 0] = 0
    # load go ids and prot ids
    goterms = []
    with open(opt.goids, "r") as go_f:
        lines = go_f.readlines()
        for line in lines:
            goterms.append(line.split("\t")[0].replace("\n", ""))
    prots = []
    with open(opt.protids, "r") as prot_f:
        lines = prot_f.readlines()
        for line in lines:
            prots.append(line.split("\t")[0])

    if opt.cutoff_hi is not None or opt.cutoff_lo is not None:

        if opt.cutoff_lo is not None:
            threshold_lo = opt.cutoff_lo
            del_grid = np.where(np.sum(annot, axis=0) < threshold_lo)[0]
            if del_grid is not None:
                annot = np.delete(annot, del_grid, axis=1)
                goterms = np.delete(goterms, del_grid)
                print("Cutting off go terms contained by less than %s proteins,"
                      " in total %s of go terms have been removed"
                      % (threshold_lo, len(del_grid)))

        if opt.cutoff_hi is not None:
            threshold_hi = opt.cutoff_hi
            del_grid = np.where(np.sum(annot, axis=0) > threshold_hi)[0]
            if del_grid is not None:
                annot = np.delete(annot, del_grid, axis=1)
                goterms = np.delete(goterms, del_grid)
                print("Cutting off go terms contained by more than %s proteins,"
                      " in total %s of go terms have been removed"
                      % (threshold_hi, len(del_grid)))


    if opt.prep:
        # go through deepNF preprocessing protocol
        Nets = loadmat(opt.input_networks)["Networks"][0]
	Nets = load_networks(Nets)

        # Compute RWR + PPMI
        for i in range(0, len(Nets)):
            print()
            print("### Computing PPMI for network:%s" % i)
            Nets[i] = RWR(Nets[i], alpha=opt.alpha)
            Nets[i] = PPMI_matrix(Nets[i])
            # print("### Writing output to file...")
            # fWrite = open('%s%s_%s_K3_alpha%d.mat' % (opt.results_path, opt.prefix, i, opt.alpha), 'wb')
            # savemat(fWrite, {'annotations': Nets[i]})
            # fWrite.close()

        savemat("%s%s_preprocessed_nets.mat" % (opt.results_path, opt.prefix),
                {"Networks":Nets}, do_compression=True)
    else:
        Nets = list(loadmat(opt.input_networks)["Networks"])
    trainer = DeepNF(annot_file=annot, networks=Nets, goids=goterms, protids=prots, batch_size=opt.batch_size,
                     cv=opt.valid_type, epochs=opt.epochs, ntrials=opt.ntrials, arch=opt.select_arch, topK=opt.topK)

    # load model if already existed
    if os.path.exists("%s%s_mda.h5" % (opt.models_path, opt.prefix)) \
            or os.path.exists("%s%s_ae.h5" % (opt.models_path, opt.prefix)):

        print("Model already Existed, Proceed to load the model")

        if os.path.exists("%s%s_mda.h5" % (opt.models_path, opt.prefix)):
            model = load_model("%s%s_mda.h5" % (opt.models_path, opt.prefix))
        else:
            model = load_model("%s%s_ae.h5" % (opt.models_path, opt.prefix))
            model = Model(inputs=model.input, outputs=model.get_layer('middle_layer').output)

    else:

        if len(Nets) > 1:
            model, history = trainer.train_mda()
            model.save(opt.models_path + opt.prefix + "_mda.h5")
        else:
            model, history = trainer.train_ae()
            model.save(opt.models_path + opt.prefix + "_ae.h5")

        plt.plot(history.history['loss'], 'o-')
        plt.plot(history.history['val_loss'], 'o-')
        plt.title('model loss')
        plt.ylabel('loss')
        plt.xlabel('epoch')
        plt.legend(['train', 'validation'], loc='upper left')
        plt.savefig(opt.models_path + opt.prefix + '_loss.png', bbox_inches='tight')

    features = model.predict(Nets)
    features = minmax_scale(features)

    scores, topKList = trainer.train_SVM(features)

    if opt.topK is not None:
        fout = open("%s%s_top%d.txt" % (opt.results_path, opt.prefix, opt.topK), "w")
        fout.write("GO ID\tProt ID\tScore\n")
        for cur_tuple in topKList:
            fout.write("%s\t%s\t%s" % (cur_tuple[0], cur_tuple[1], cur_tuple[2]))
        fout.close()

        # measures
    measures = ['m-aupr_avg', 'm-aupr_std', 'M-aupr_avg', 'M-aupr_std',
                'F1_avg', 'F1_std', 'acc_avg', 'acc_std']

    fout = open("%s%s_performance_score.txt" % (opt.results_path, opt.prefix), "w")
    fout.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ('m-aupr_avg', 'm-aupr_std', 'M-aupr_avg', 'M-aupr_std',
                'F1_avg', 'F1_std', 'acc_avg', 'acc_std'))
    for m in measures:
        fout.write('%0.5f ' % (scores[m]))
    fout.write('\n')
    fout.close()


if __name__ == "__main__":
    run()
