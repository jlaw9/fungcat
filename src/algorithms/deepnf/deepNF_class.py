import os
os.environ["KERAS_BACKEND"] = "tensorflow"
from keras.models import Model, load_model
from validation import cross_validation, temporal_holdout
from deepNF import build_MDA, build_AE
from keras.callbacks import EarlyStopping
import matplotlib as mpl
mpl.use('Agg')
import numpy as np


class DeepNF:

    def __init__(self, networks, annot_file, goids, protids, epochs, ntrials, batch_size,
                 topK, cv, arch):

        # all possible combinations for architectures
        if arch is None:
            self.arch = [6*2500, 1200, 6*2500]
        else:
            self.arch = arch


        nets = []
        input_dims = []

        for net in networks:
            nets.append(net)
            input_dims.append(net.shape[1])
        self.nets = nets
        self.input_dims = input_dims
        self.goterms = goids
        self.prots = protids
        self.ntrials = ntrials
        self.epochs = epochs
        self.batch_size = batch_size
        self.cv = cv
        self.annot = annot_file
        self.topK = topK

    def train_mda(self):

        print "### [Model] Running for architecture: ", self.arch

        model = build_MDA(self.input_dims, self.arch)
        history = model.fit(self.nets, self.nets, epochs=self.epochs, batch_size=self.batch_size, shuffle=True,
                            validation_split=0.1,
                            callbacks=[EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=2)])
        mid_model = Model(inputs=model.input,
                          outputs=model.get_layer('middle_layer').output)

        return mid_model, history

    def train_ae(self):

        print "### [Model] Running for architecture: ", self.arch
        print "### [Model 1] Running for network1 "
        model = build_AE(self.input_dims[0], self.arch)
        history = model.fit(self.nets, self.nets, epochs=self.epochs, batch_size=self.batch_size, shuffle=True,
                            validation_split=0.1,
                            callbacks=[EarlyStopping(monitor='val_loss', min_delta=0.0001, patience=2)])
        return model, history

    def train_SVM(self, features):

        print "### Running for SVM prediction:"
        if self.cv:
            perf = cross_validation(features, self.annot,
                                    n_trials=self.ntrials, goterms=self.goterms)
            topK_list = []
            if self.topK is not None:

                if self.topK > perf['pr_goterms'].shape[1]:
                    print "k of %s is bigger than the given size of protein: %s" % \
                          (self.topK, perf['pr_goterms'].shape[1])
                else:
                    # topK = np.sort(np.delete(sorted, np.s_[k::], axis=1), axis=1)
                    pair_matrix = perf['pr_goterms'].T
                    row_num, column_num = pair_matrix.shape
                    for row in np.arange(row_num):
                        indices = np.arange(column_num)
                        # sort two lists together
                        cur_goterm, indices = zip(*sorted(zip(pair_matrix[row], indices), reverse=True))
                        cur_goterm = np.delete(cur_goterm, np.s_[self.topK::])
                        indices = np.delete(indices, np.s_[self.topK::])
                        goterm_name = self.goterms[row]
                        for column in np.arange(self.topK):
                            topK_list.append((goterm_name, self.prots[indices[column]], cur_goterm[column]))

            return perf, topK_list
        return
