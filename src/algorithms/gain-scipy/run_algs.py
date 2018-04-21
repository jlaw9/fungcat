
# Quick script to run/test the algorithms
print("Importing libraries")

from optparse import OptionParser
from collections import defaultdict
import os
import sys
from tqdm import tqdm
import itertools
import utils.file_utils as utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fungcat_settings as f_settings
import alg_utils
import sinksource
#import sinksource_topk_ub
import sinksource_ripple
import sinksource_squeeze
#import pandas as pd
import networkx as nx
import scipy
import numpy as np


class Alg_Runner:
    """
    Base class for running algorithms
    """
    def __init__(self, versions, exp_name, pos_neg_files, goterms=None, 
            algorithms=["sinksource"], eps_list=[0.0001], k_list=[100], t_list=[2], s_list=[50], a_list=[0.8], 
            deltaUBLB_list=[None], epsUB_list=[0], taxon=None, num_pred_to_write=100, cross_validation_folds=None, 
            forcealg=False, verbose=False):
        """
        *eps*: Convergence cutoff for sinksource and sinksourceplus
        """
        self.versions = versions
        self.exp_name = exp_name
        self.goterms = goterms
        self.pos_neg_files = pos_neg_files
        self.algorithms = algorithms
        self.eps_list = eps_list
        self.k_list = k_list
        self.t_list = t_list
        self.s_list = s_list
        self.a_list = a_list
        self.deltaUBLB_list = deltaUBLB_list
        self.epsUB_list = epsUB_list
        self.taxon = taxon
        self.cross_validation_folds = cross_validation_folds
        self.num_pred_to_write = num_pred_to_write
        self.forcealg = forcealg
        self.verbose = verbose

        # will be instantiated later
        self.int2node = {}
        self.node2int = {}

    def parse_pos_neg_matrix(self, pos_neg_file):
        print("Reading positive and negative annotations for each protein from %s" % (pos_neg_file))
        goid_prots = defaultdict(set)
        goid_neg = defaultdict(set)
        num_prots = 0
        # TODO write the BP and MF terms to the same file
        # TODO possibly use pickle
        goid_prots_file = pos_neg_file.replace(".tsv", "-pos-list.txt")
        goid_neg_file = pos_neg_file.replace(".tsv", "-neg-list.txt")
        if os.path.isfile(goid_prots_file):
            goid_prots = {goid: set(prots.split(',')) for goid, prots in utils.readColumns(goid_prots_file, 1, 2)}
            goid_neg = {goid: set(prots.split(',')) for goid, prots in utils.readColumns(goid_neg_file, 1, 2)}
        else:
            #pos_neg = df.read_csv(pos_neg_file)
            # for each GO term, get the set of positives and negatives
            with open(pos_neg_file, 'r') as f:
                goterms = f.readline().rstrip().split('\t')[1:]
                for line in f:
                    line = line.rstrip().split('\t')
                    num_prots += 1
                    prot = line[0]
                    vals = line[1:]
                    for i in range(0, len(vals)):
                        val = int(vals[i])
                        if val == 1:
                            goid_prots[goterms[i]].add(prot)
                        elif val == -1:
                            goid_neg[goterms[i]].add(prot)

            print("Writing %s for easier access" % (goid_prots_file))
            with open(goid_prots_file, 'w') as out:
                out.write(''.join(["%s\t%s\n" % (goid, ','.join(prots)) for goid, prots in goid_prots.items()]))

            with open(goid_neg_file, 'w') as out:
                out.write(''.join(["%s\t%s\n" % (goid, ','.join(prots)) for goid, prots in goid_neg.items()]))

        print("\t%d GO terms, %d prots" % (len(goid_prots), num_prots))

        return goid_prots, goid_neg


    def parse_pos_neg_file(self, pos_neg_file):
        print("Reading positive and negative annotations for each protein from %s" % (pos_neg_file))
        goid_prots = {}
        goid_neg = {}
        all_prots = set()
        # TODO possibly use pickle
        if os.path.isfile(pos_neg_file):
            for goid, pos_neg_assignment, prots in utils.readColumns(pos_neg_file, 1,2,3):
                prots = set(prots.split(','))
                if int(pos_neg_assignment) == 1:
                    goid_prots[goid] = prots
                elif int(pos_neg_assignment) == -1:
                    goid_neg[goid] = prots

                all_prots.update(prots)

        print("\t%d GO terms, %d prots" % (len(goid_prots), len(all_prots)))

        return goid_prots, goid_neg


    def setup_network(self, network_file):
        # TODO update reading the initial network file
        #normalized_net_file = network_file.replace(".txt", "-normalized.txt")
        sparse_net_file = network_file.replace('.txt', '-normalized.npz')
        node2int_file = network_file.replace(".txt", "-node2int.txt")
        if os.path.isfile(sparse_net_file):
            print("Reading network from %s" % (sparse_net_file))
            P = scipy.sparse.load_npz(sparse_net_file)
#        else:
#            print("Reading network from %s" % (normalized_net_file))
#            # takes more memory and time to read the file into memory and then loop through it
#            #lines = utils.readColumns(normalized_net_file, 1, 2, 3)
#            #H.add_weighted_edges_from([(int(u),int(v),float(w)) for u,v,w in lines])
#            # TODO store the network as a 
#            i, j, val = np.loadtxt(normalized_net_file).T
#            P = scipy.sparse.coo_matrix((val, (i, j))).tocsr()
#            scipy.sparse.save_npz(sparse_net_file, P)
            #print("\t%d nodes and %d edges" % (H.number_of_nodes(), H.number_of_edges()))
            print("Reading node names from %s" % (node2int_file))
            node2int = {n: int(n2) for n, n2 in utils.readColumns(node2int_file, 1, 2)}
            int2node = {int(n): n2 for n, n2 in utils.readColumns(node2int_file, 2, 1)}

        else:
            print("Reading network from %s" % (network_file))
            # TODO read network and convert to sparse matrix without networkx
    #            lines = utils.readColumns(network_file, 1, 2, 3)
    #            nodes = set(n for u,v,w in lines for n in (u,v))
    #            # map the nodes to integers
    #            alg_utils.convert_nodes_to_int(nodes)
            G = nx.Graph()
            #G.add_weighted_edges_from([(u,v,float(w)) for u,v,w in lines])
            #edges = {}
            #nodes = set()
            with open(network_file, 'r') as f:
                for line in f:
                    if line[0] == "#":
                        continue
                    u,v,w = line.rstrip().split('\t')[:3]
                    G.add_edge(u,v,weight=float(w))
                    #edges[(u,v)] = float(w)
                    #nodes.update(set([u,v]))
            print("\t%d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))

            print("\trelabeling node IDs with integers")
            G, node2int, int2node = alg_utils.convert_nodes_to_int(G)
            print("\twriting node2int labels to %s" % (node2int_file))
            with open(node2int_file, 'w') as out:
                out.write(''.join(["%s\t%s\n" % (int2node[n], n) for n in sorted(int2node)]))

            print("\tconverting to a scipy sparse matrix")
            W = nx.to_scipy_sparse_matrix(G, nodelist=sorted(G.nodes()))

            print("\tnormalizing edge weights by dividng each edge weight by tail node's out_degree")
            P = alg_utils.normalizeGraphEdgeWeights(W)

            print("\twriting sparse normalized network to %s" % (sparse_net_file))
            scipy.sparse.save_npz(sparse_net_file, P)

            #print("\twriting normalized network to %s" % (normalized_net_file))
            #with open(normalized_net_file, 'w') as out:
            #    out.write(''.join(["%s\t%s\t%s\n" % (u,v,str(data['weight'])) for u,v,data in H.edges(data=True)]))
        print("\t%d nodes and %d edges" % (P.shape[0], len(P.data)))
        # as a sanity check print the first node's edges
        #print(P[0])

        return P, node2int, int2node


    def main(self):
        version_params_results = {}

        #version = "2017_10-seq-sim"
        #versions = ["2017_10-string"]
        for version in self.versions:
            INPUTSPREFIX, RESULTSPREFIX, network_file, selected_strains = f_settings.set_version(version)
            self.INPUTSPREFIX = INPUTSPREFIX
            self.RESULTSPREFIX = RESULTSPREFIX
            #gain_file = f_settings.GOA_ALL_FUN_FILE
            # start with a single species
            #self.taxon = "208964"
            if self.taxon is not None:
                network_file = f_settings.STRING_TAXON_UNIPROT % (self.taxon, self.taxon, f_settings.STRING_CUTOFF)
                gain_file = f_settings.FUN_FILE % (self.taxon, self.taxon)

            P, node2int, int2node = self.setup_network(network_file)
            self.node2int = node2int
            self.int2node = int2node

            # get the positives and negatives from the matrix
            goid_pos_neg = {}
            all_goid_prots = {}
            all_goid_neg = {}
            for pos_neg_file in self.pos_neg_files:
                #goid_prots, goid_neg = self.parse_pos_neg_matrix(self.pos_neg_file)
                goid_prots, goid_neg = self.parse_pos_neg_file(pos_neg_file)
                all_goid_prots.update(goid_prots)
                all_goid_neg.update(goid_neg)

            print(len(self.goterms))
            if self.goterms is None:
                self.goterms = set(goid_pos_neg.values())

            for goterm in sorted(self.goterms):
                if goterm not in all_goid_prots:
                    print("Skipping %s. Not in pos-neg-file" % (goterm))
                    continue
                # covert the positives to integers to match the graph
                positives = np.asarray([node2int[n] for n in all_goid_prots[goterm] if n in node2int])
                negatives = np.asarray([node2int[n] for n in all_goid_neg[goterm] if n in node2int])
                goid_pos_neg[goterm] = {'pos':positives, 'neg':negatives}

            #del all_goid_prots
            #del all_goid_neg

            # organize the results by algorithm, then algorithm parameters, then GO term
            params_results = {}
            for alg in self.algorithms:
                out_dir = "%s/all/%s/%s" % (self.RESULTSPREFIX, alg, self.exp_name)
                utils.checkDir(out_dir)

                # instead of storing the results for each GO term in a separate file,
                # store the predictions for all GO terms in the same file for each algorithm and parameter set
                # For now, each run is appending to the file, rather than writing everything at the end
                # TODO check if the file already exists. 
                # If it does, load it and check if the results for the GO term already exist, 
                # and proceed accordingly
                #out_file = "%s/%s"

                #out_pref = "%s/%s-" % (out_dir, goterm.replace(':',''))
                if self.num_pred_to_write == 0:
                    out_pref = None
                else:
                    out_pref = "%s/all-pred-" % (out_dir)

                # now run the algorithm with all combinations of parameters
#                    if alg in ["sinksource", "sinksourceplus"]:
                curr_params_results = self.run_alg_with_params(
                        P, alg, goid_pos_neg, out_pref=out_pref)

#                    elif alg in ["sinksource-ripple", "sinksourceplus-ripple"]:
#                        scores, curr_params_results = run_topk_with_params(
#                                P, alg, goterm, out_pref, positives, negatives,
#                                k_list, t_list, s_list, a_list, deltaUBLB_list, 
#                                epsUB_list, forced=forcealg)
                
                params_results.update(curr_params_results)

                if self.cross_validation_folds is not None:
                    cv_params_results = run_cross_validation(H, alg, goterm, positives, negatives, out_pref, folds=cross_validation_folds, a=0.8, eps=0.0001)
                    params_results.update(cv_params_results)

                if self.verbose:
                    print(params_results_to_table(params_results))

            version_params_results[version] = params_results

        if self.verbose:
            for version in self.versions:
                print(version)
                table = params_results_to_table(params_results)
                out_dir = "outputs/viz/%s" % (self.exp_name)
                utils.checkDir(out_dir)
                out_file = "%s/params-results.tsv" % (out_dir)
                print("writing %s" % (out_file))
                with open(out_file, 'w') as out:
                    out.write(table)
                print(table)

        return

    ## TODO test the difference between the GAIN pos/neg and our method of assigning positives and negatives


    def run_alg_with_params(self, P, alg, goid_pos_neg,
                            out_pref=None):
        """
        """
        params_results = {}

        # TODO run these separately because they have different parameters
        if alg == "sinksourceplus" or alg == "sinksource":
            params_results = self.run_ss_with_params(
                    P, alg, goid_pos_neg,
                    out_pref=out_pref)
        elif alg in ["sinksource-ripple", "sinksourceplus-ripple", 
                'sinksource-squeeze', 'sinksourceplus-squeeze']:
            params_results = self.run_topk_with_params(
                    P, alg, goid_pos_neg,
                    out_pref=out_pref)
                    #out_pref=out_pref, forced=forced)

        return params_results


    def run_ss_with_params(self, P, alg, goid_pos_neg,
                           out_pref=None):
        all_params_results = {} 
        #for a, eps in tqdm(list(itertools.product(*[a_list, eps_list]))):
        for a, eps in itertools.product(*[self.a_list, self.eps_list]):
            out_file = None
            if out_pref is not None:
                out_file = "%sa%s-eps%s.txt" % (
                        out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'))
                #header = "#%s\tGO term: %s\t%d positives\t%d negatives\ta=%s\teps=%s\n" \
                #        % (alg, goterm, len(positives), len(negatives), str(a), str(eps))

            params_results = self.run_alg_on_goterms(P, alg, goid_pos_neg, 
                    out_file=out_file, a=a, eps=eps)
            all_params_results.update(params_results)

            if self.verbose:
                print(params_results_to_table(params_results))
        return all_params_results


    def run_topk_with_params(self, P, alg, goid_pos_neg,
                             out_pref=None):
        all_params_results = {}

        # Run using all combinations of parameters
        # s and t won't change the results, just the amount of time, so loop through those separately
        total_params_list = list(itertools.product(*[self.k_list, self.t_list, self.s_list, self.a_list, self.epsUB_list]))
        #params_list = list(itertools.product(*[k_list, a_list, deltaUBLB_list]))
        #t_s_list = list(itertools.product(*[t_list, s_list]))
        print("Running %d combinations of parameters" % (len(total_params_list)))
        for k, a, epsUB in itertools.product(*[self.k_list, self.a_list, self.epsUB_list]):

            if "ripple" in alg:
                for t, s in itertools.product(*[self.t_list, self.s_list]):
                    out_file = None
                    if out_pref is not None:
                        out_file = "%sk%d-t%d-s%d-a%s-eps%s-scipy.txt" % (out_pref, k, t, s, 
                                str(a).replace('.', '_'), str(epsUB).replace('.', '_'))

                    #print("Running %s with %d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, deltaUBLB=%s for GO term %s" \
                    #        % (alg, len(positives), len(negatives), k, t, s, str(a), str(deltaUBLB), goterm))

                    params_results = self.run_alg_on_goterms(P, alg, goid_pos_neg, out_file=out_file,
                            a=a, k=k, t=t, s=s, epsUB=epsUB)
                    all_params_results.update(params_results)
                        # also keep track of the time it takes for each of the parameter sets

            elif 'squeeze' in alg:
                out_file = None
                if out_pref is not None:
                    out_file = "%sk%d-a%s-epsUB%s.txt" % (out_pref, k, str(a).replace('.', '_'),
                                                            str(epsUB).replace('.', '_'))
                params_results = self.run_alg_on_goterms(P, alg, goid_pos_neg, out_file=out_file,
                        a=a, k=k, epsUB=epsUB)
                all_params_results.update(params_results)

            if self.verbose:
                print(params_results_to_table(params_results))

        return all_params_results


    def run_alg_on_goterms(self, P, alg, goid_pos_neg, out_file=None, 
            a=0.8, eps='-', k='-', t='-', s='-', epsUB=0):
                # saves a lot of time keeping the file handle open
        params_results = {}
        try:
            if out_file is not None:
                if self.forcealg is False and os.path.isfile(out_file):
                    print("%s already exists. Use --forcealg to overwrite" % (out_file))
                    return params_results

                file_handle = open(out_file, 'w')
                file_handle.write("#goterm\tprot\tscore\n")
                #file_handle.write("#%s\tk=%d\tt=%d\ts=%d\ta=%s\tepsUB=%s\n" \
                #        % (alg, k, t, s, str(a), str(epsUB)))

            print("Running %s for %d goterms. Writing to %s" % (alg, len(goid_pos_neg), out_file))
            for goterm in tqdm(goid_pos_neg):
                positives, negatives = goid_pos_neg[goterm]['pos'], goid_pos_neg[goterm]['neg']
                num_unk = P.shape[0] - len(positives)
                if alg in ['sinksource', 'sinksource-squeeze', 'sinksource-ripple']:
                    num_unk -= len(negatives) 
                len_N = P.shape[0]

                # TODO streamline calling the correct function. They all take the same parameters
                # This uses the same UB as Ripple
                if alg == "sinksourceplus":
                    scores, time, iters, comp = sinksource.runSinkSource(
                            P, positives, negatives=None, max_iters=1000, delta=eps, a=a)
                elif alg == "sinksource":
                    scores, time, iters, comp = sinksource.runSinkSource(
                            P, positives, negatives=negatives, max_iters=1000, delta=eps, a=a)
                elif alg == "sinksourceplus-squeeze":
                    R, scores, time, iters, comp, max_ds = sinksource_squeeze.runSinkSourceSqueeze(
                            P, positives, k=k, a=a, epsUB=epsUB, verbose=self.verbose)
                elif alg == "sinksource-squeeze":
                    R, scores, time, iters, comp, max_ds = sinksource_squeeze.runSinkSourceSqueeze(
                            P, positives, negatives=negatives, k=k, a=a, epsUB=epsUB, verbose=self.verbose)
                elif alg == "sinksourceplus-ripple":
                    R, scores, time, iters, comp, len_N, max_ds = sinksource_ripple.runSinkSourceRipple(
                            P, positives, k=k, t=t, s=s, a=a, epsUB=epsUB)
                # This uses the same UB as Ripple, but with negatives
                elif alg == "sinksource-ripple":
                    R, scores, time, iters, comp, len_N, max_ds = sinksource_ripple.runSinkSourceRipple(
                            P, positives, negatives, k=k, t=t, s=s, a=a, epsUB=epsUB)

                tqdm.write("\t%s converged after %d iterations " % (alg, iters) + \
                        "(%0.2f sec) for goterm %s" % (time, goterm))

                # also keep track of the time it takes for each of the parameter sets
                params_key = (alg, goterm, len(positives), num_unk, a, eps, k, t, s, epsUB)
                params_results[params_key] = (time, iters, comp, len_N)

                if out_file is not None:
                    # convert the nodes back to their names, and make a dictionary out of the scores
                    scores = {self.int2node[n]:scores[n] for n in scores}
                    self.write_scores_to_file(scores, goterm=goterm, file_handle=file_handle,
                            num_pred_to_write=self.num_pred_to_write)

                # TODO move this somewhere else
                # plot the max_ds
                self.plot_max_ds(max_ds)
        except:
            if out_file is not None:
                file_handle.close()
            raise

        if out_file is not None:
            file_handle.close()
            print("Finished running %s for %d goterms. Wrote to %s" % (alg, len(goid_pos_neg), out_file))

        return params_results


    def plot_max_ds(self, max_d_list):
#import matplotlib.pyplot as plt
#plt.plot(max_d_list)
#plt.yaxis("Maximum score change")
#plt.xaxis("iteration")
        out_file = "outputs/viz/max-d.txt"
        print("writing %s" % out_file)
        with open(out_file, 'w') as out:
            out.write('\n'.join(str(x) for x in max_d_list))
        #plt.savefig(out_file)


    def write_scores_to_file(self, scores, goterm='', out_file=None, file_handle=None,
            num_pred_to_write=100, header="", append=True):
        """
        *scores*: dictionary of node_name: score
        *num_pred_to_write*: number of predictions (node scores) to write to the file 
            (sorted by score in decreasing order). If -1, all will be written
        """

        if num_pred_to_write == -1:
            num_pred_to_write = len(scores) 

        if out_file is not None:
            if append:
                print("Appending %d scores for goterm %s to %s" % (num_pred_to_write, goterm, out_file))
                out_type = 'a'
            else:
                print("Writing %d scores for goterm %s to %s" % (out_file))
                out_type = 'w'

            file_handle = open(out_file, out_type)
        elif file_handle is None:
            print("Warning: out_file and file_handle are None. Not writing scores to a file")
            return 

        # write the scores to a file
        # for now, write the top k node's results
        # TODO Possibly write the results to a JSON file so it's easier to load
        file_handle.write(header)
        for n in sorted(scores, key=scores.get, reverse=True)[:num_pred_to_write]:
            file_handle.write("%s\t%s\t%s\n" % (goterm, n, str(scores[n])))
        return


def run_cross_validation(H, alg, goterm, positives, negatives, out_pref, folds=5, a=0.8, eps=0.0001):
    from sklearn.model_selection import KFold

    # compare how long it takes to run each from scratch vs using the previous run's scores
    global int2node
    global node2int
    params_results = {}
    overall_time_normal = 0
    overall_time_using_predictions = 0

    if alg != "sinksourceplus":
        print("Not yet implemented")
        return
    if alg == "sinksourceplus":
        # this should already have been created from the regular run
        out_file = "%sa%s-eps%s.txt" % (out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'))
        print("\tReading scores from %s" % (out_file))
        prediction_scores = {node2int[goterm]: float(score) for goterm, score in utils.readColumns(out_file, 1, 2)}
        #predition_scores, time, iters = sinksource.runSinkSource(H, positives, negatives=None, max_iters=1000, delta=eps, a=a)

    # make the positives and negatives into numpy arrays
    positives = np.array(list(positives))
    #negatives = np.array(negatives)
    kf = KFold(n_splits=folds, shuffle=True)
    kf.get_n_splits(positives)
    fold = 0
    orig_num_iters = []
    scratch_num_iters = []
    for train_idx, test_idx in kf.split(positives):
        fold += 1
        pos_train, pos_test = set(positives[train_idx]), set(positives[test_idx])
        print(pos_train)
        print("Current fold: %d train positives, %d test" % (len(pos_train), len(pos_test)))
        if alg == "sinksourceplus":
            scores, time, iters = sinksource.runSinkSource(H, pos_train, negatives=None, max_iters=1000, delta=eps, a=a, scores=prediction_scores.copy())
            orig_num_iters.append(iters)
            overall_time_using_predictions += time
            print("Time (min) to run from original scores: %0.2f" % (time / float(60)))
#            out_file = "%sa%s-eps%s-fold%d-orig.txt" % (out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'), fold)
#            print("Writing %s" % (out_file))
#            # write the scores to a file
#            with open(out_file, 'w') as f:
#                for n in sorted(scores, key=scores.get, reverse=True):
#                    f.write("%s\t%0.6f\n" % (int2node[n], scores[n]))

            scores, time, iters = sinksource.runSinkSource(H, pos_train, negatives=None, max_iters=1000, delta=eps, a=a)
            scratch_num_iters.append(iters)
            overall_time_normal += time
            print("Time (min) to run from scratch: %0.2f" % (time / float(60)))
#            out_file = "%sa%s-eps%s-fold%d-scratch.txt" % (out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'), fold)
#            print("Writing %s" % (out_file))
#            # write the scores to a file
#            with open(out_file, 'w') as f:
#                for n in sorted(scores, key=scores.get, reverse=True):
#                    f.write("%s\t%0.6f\n" % (int2node[n], scores[n]))
            # TODO evaluate the cross-validation
    print("Total time (min) to run from scratch: %0.2f" % (overall_time_normal / float(60)))
    print("Total time (min) to run from original scores: %0.2f" % (overall_time_using_predictions / float(60)))
    print("Iterations from scratch: %s" % (','.join([str(x) for x in scratch_num_iters])))
    print("Iterations from original scores: %s" % (','.join([str(x) for x in orig_num_iters])))
    time_diff = overall_time_normal - overall_time_using_predictions
    avg_iter_diff = np.mean(scratch_num_iters) - np.mean(orig_num_iters)
    params_key = (alg, goterm, 'folds=5', '-', '-', a, 'eps=%s'%str(eps))
    params_results[params_key] = (time_diff, avg_iter_diff, '-')

    return params_results

#    # TODO run these separately because they have different parameters
#    if alg == "sinksourceplus" or alg == "sinksource":
#
#        if alg == "sinksourceplus":
#            scores, time, iters = sinksource.runSinkSource(H, pos_train, negatives=None, max_iters=1000, delta=eps, a=a)
#        if alg == "sinksource":
#            scores, time, iters = sinksource.runSinkSource(H, pos_train, negatives=negatives, max_iters=1000, delta=eps, a=a)
##        params_key = (alg, goterm, '-', '-', '-', a, 'eps=%s'%str(eps))
##        params_results[params_key] = (time, iters, '-')
##
##        print(params_results_to_table(params_results))
##
##        return params_results
#
#    else:
#        print("Not yet implemented")
#        return


def params_results_to_table(params_results):
    # print out a table of the results for each parameter set
    # print out the table after each iteration so the results are still viewable if the script quits earl or something
    results = "\t".join(['alg', 'goterm', '# pos', '# unk', 'a', 'k', 't', 's', 'eps', 'epsUB', 'time', 'iters', '# comp', 'len_N']) + '\n'
    for params_key in sorted(params_results):
        alg, goterm, num_pos, num_unk, a, eps, k, t, s, epsUB = params_key
        time, iters, total_comp, len_N = params_results[params_key]
        results += "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%0.2f\t%d\t%0.2e\t%s\n" % \
                (alg, goterm, num_pos, num_unk, str(a), str(k), str(t), str(s), str(eps), str(epsUB), 
                        time, iters, total_comp, str(len_N))
    return results


ALGORITHMS = [
    "sinksourceplus-ripple",  # This uses the same UB as Ripple
    "sinksource-ripple",  # This uses the same UB as Ripple, but with negatives
#    "sinksource-topk-ub",  # This computes UB and LB using iteration
#    "sinksourceplus-topk-ub",  # This uses the same UB as sinksource-ub-topk
    "sinksourceplus-squeeze",
    "sinksource-squeeze",
    "sinksourceplus",
    "sinksource",
    ]

def parse_args(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--version',type='string', action='append',
                      help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    parser.add_option('-A', '--algorithm', action="append",
                      help="Algorithm for which to get predictions. Default is all of them. Options: '%s'" % ("', '".join(ALGORITHMS)))
    parser.add_option('', '--exp-name', type='string',
                      help="Experiment name to use when running GAIN.")
    parser.add_option('', '--pos-neg-file', type='string', action='append',
                      help="File containing positive and negative examples for each GO term")
    parser.add_option('', '--only-functions', type='string',
                      help="Run GAIN using only the functions in a specified file (should be in XX format i.e., without the leading GO:00).")
    parser.add_option('-G', '--goterm', type='string', action="append",
                      help="Specify the GO terms to use (should be in GO:00XX format)")
    parser.add_option('-k', '--k', type=int, action="append",
                      help="Top-k for Ripple, Default=200")
    parser.add_option('-t', '--t', type=int, action="append",
                      help="t parameter for Ripple. Default=2")
    parser.add_option('-s', '--s', type=int, action="append",
                      help="s parameter for Ripple. Default=200")
    parser.add_option('-a', '--alpha', type=float, action="append",
                      help="Alpha insulation parameter. Default=0.8")
    parser.add_option('', '--eps', type=float, action="append",
                      help="Stopping criteria for SinkSource")
    parser.add_option('-e', '--epsUB', type=float, action="append",
                      help="Parameter to return the top-k if all other nodes have an UB - epsUB < the kth node's LB. Default=0")
    parser.add_option('-W', '--num-pred-to-write', type='int', default=100,
                      help="Number of predictions to write to the file. If 0, none will be written. If -1, all will be written. Default=100")
    parser.add_option('-C', '--cross-validation-folds', type='int',
                      help="Perform cross validation using the specified # of folds. Usually 5")
    parser.add_option('', '--forcealg', action="store_true", default=False,
                      help="Force re-running algorithms if the output files already exist")
    parser.add_option('', '--verbose', action="store_true", default=False,
                      help="Print additional info about running times and such")

    #parser.add_option('', '--goid', type='string', metavar='STR',
    #                  help='GO-term ID for which annotations and precitions will be posted')

    (opts, args) = parser.parse_args(args)

    if opts.exp_name is None or opts.pos_neg_file is None:
        print("--exp-name, --pos-neg-file, required")
        sys.exit(1)

    # if neither are provided, just use all GO terms in the pos/neg file
#    if opts.goterm is None and opts.only_functions is None:
#        print("--goterm or --only_functions required")
#        sys.exit(1)

    if opts.version is None:
        opts.version = ["2017_10-seq-sim"]

    if opts.algorithm is None:
        opts.algorithm = ALGORITHMS

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)

    #if opts.algorithm not in f_settings.ALGORITHM_OPTIONS:
    #    print "--algorithm %s is not a valid algorithm name" % (opts.algorithm)
    #    sys.exit(1)

    # TODO
    for alg in opts.algorithm:
        if alg not in ALGORITHMS:
            print("ERROR: '%s' not a valid algorithm name. Algorithms are: '%s'." % (alg, ', '.join(ALGORITHMS)))
            sys.exit(1)

    opts.k = opts.k if opts.k is not None else [200]
    opts.t = opts.t if opts.t is not None else [2]
    opts.s = opts.s if opts.s is not None else [200]
    opts.alpha = opts.alpha if opts.alpha is not None else [0.8]
    # default for deltaUBLB is None
    opts.eps = opts.eps if opts.eps is not None else [0.0001]
    opts.epsUB = opts.epsUB if opts.epsUB is not None else [0]

    return opts


if __name__ == "__main__":
    #versions = ["2017_10-seq-sim", "2017_10-seq-sim-x5-string"]
    opts = parse_args(sys.argv)

    goterms = None
    if opts.only_functions is not None:
        only_functions = utils.readItemSet(opts.only_functions, 1)
        goterms = set(["GO:" + "0"*(7-len(str(x))) + str(x) for x in only_functions])
    if opts.goterm is not None:
        goterms = set() if goterms is None else goterms
        goterms.update(set(opts.goterm))

    alg_runner = Alg_Runner(opts.version, opts.exp_name,
         opts.pos_neg_file, goterms, opts.algorithm,
         k_list=opts.k, t_list=opts.t, s_list=opts.s, a_list=opts.alpha,
         eps_list=opts.eps, epsUB_list=opts.epsUB,
         num_pred_to_write=opts.num_pred_to_write,
         cross_validation_folds=opts.cross_validation_folds,
         forcealg=opts.forcealg, verbose=opts.verbose)
    alg_runner.main()
    #main()
