
# Quick script to run/test the algorithms
print("Importing libraries")

from optparse import OptionParser
import os
import sys
from tqdm import tqdm
import utils.file_utils as utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fungcat_settings as f_settings
import networkx as nx
import alg_utils
import sinksource
import sinksourceplus_topk
import sinksource_ripple
from collections import defaultdict
#import pandas as pd


ALGORITHMS = [
    "sinksourceplus-ripple",  # This uses the same UB as Ripple
    "sinksource-ripple",  # This uses the same UB as Ripple, but with negatives
    "sinksource-topk-ub",  # This computes UB and LB using iteration
    "sinksourceplus-topk-ub",  # This uses the same UB as sinksource-ub-topk
    ]


def parse_gain_file(gain_file, goterms=None):
    print("Reading annotations from GAIN file %s. Assuming annotations have already been propogated up the GO DAG" % (gain_file))
    goid_prots = defaultdict(set)

    with open(gain_file, 'r') as f:
        # columns: 'orf', 'goid', 'hierarchy', 'evidencecode', 'annotation type' (described here: http://bioinformatics.cs.vt.edu/~murali/software/biorithm/gain.html)
        header_line = f.readline().rstrip().split('\t')
        # *orf* A systematic name for the gene.
        orf_index = header_line.index('orf')
        # *goid* The ID of the GO function. You can leave in the "GO:0+" prefix for a function. GAIN will strip it out.
        goid_index = header_line.index('goid')
        goname_index = header_line.index('goname')
        # *hierarchy* The GO category the function belongs to ("c", "f", and "p")
        hierarchy_index = header_line.index('hierarchy')
        # *evidencecode* The evidence code for an annotation
        evidencecode_index = header_line.index('evidencecode')
        # *annotation type* The value is either 1 (annotated), 0 (unknown), or -1 (not annotated)
        ann_index = header_line.index('annotation type')

        for line in f:
            line = line.rstrip().split('\t')
            prot = line[orf_index]
            goid = line[goid_index]
            goname = line[goname_index]
            hierarchy = line[hierarchy_index]
            evidencecode = line[evidencecode_index]
            annotation_type = line[ann_index]

            # only track the annotations that are positive
            if annotation_type != "1":
                continue

            #hierarchies[goid] = hierarchy
            #gonames[goid] = goname
            if goterms is None or goid in goterms:
                goid_prots[goid].add(prot)

    return goid_prots


def parse_pos_neg_matrix(pos_neg_file):
    print("Reading positive and negative annotations for each protein from %s" % (pos_neg_file))
    goid_prots = defaultdict(set)
    goid_neg = defaultdict(set)
    num_prots = 0
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

    print("\t%d GO terms, %d prots" % (len(goid_prots), num_prots))

    return goid_prots, goid_neg


#def main(version, pos_neg_file, algorithms=["sinksource-ripple"], taxon=None, k_list=[100], t_list=[2], s_list=[50], a_list=[0.8], eps_list=[None]):
def main():
    version = "2017_10-seq-sim"
    #version = "2017_10-string"
    INPUTSPREFIX, RESULTSPREFIX, network_file, selected_strains = f_settings.set_version(version)
    gain_file = f_settings.GOA_ALL_FUN_FILE
    # start with a single species
    taxon = None
    #taxon = "208964"
    if taxon is not None:
        network_file = f_settings.STRING_TAXON_UNIPROT % (taxon, taxon, f_settings.STRING_CUTOFF)
        gain_file = f_settings.FUN_FILE % (taxon, taxon)

    # TODO write a normalized version of the graph so I don't have to repeat this step each time
    normalized_net_file = network_file.replace(".txt", "-normalized.txt")
    node2int_file = network_file.replace(".txt", "-node2int.txt")
    if os.path.isfile(normalized_net_file):
        H = nx.DiGraph()
        print("Reading network from %s" % (normalized_net_file))
        lines = utils.readColumns(normalized_net_file, 1, 2, 3)
        H.add_weighted_edges_from([(int(u),int(v),float(w)) for u,v,w in lines])
        print("\t%d nodes and %d edges" % (H.number_of_nodes(), H.number_of_edges()))
        print("Reading node names from %s" % (node2int_file))
        node2int = {n: int(n2) for n, n2 in utils.readColumns(node2int_file, 1, 2)}
        int2node = {int(n): n2 for n, n2 in utils.readColumns(node2int_file, 2, 1)}

    else:
        G = nx.Graph()
        print("Reading network from %s" % (network_file))
        lines = utils.readColumns(network_file, 1, 2, 3)
        G.add_weighted_edges_from([(u,v,float(w)) for u,v,w in lines])
        print("\t%d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))

        print("\trelabeling node IDs with integers")
        H, node2int, int2node = alg_utils.convert_labels_to_int(G)

        print("\tnormalizing edge weights by each node's out_degree")
        H = alg_utils.normalizeGraphEdgeWeights(H)

        print("\twriting normalized network to %s" % (normalized_net_file))
        with open(normalized_net_file, 'w') as out:
            out.write(''.join(["%s\t%s\t%s\n" % (u,v,str(data['weight'])) for u,v,data in H.edges(data=True)]))
        print("\twriting node2int labels to %s" % (node2int_file))
        with open(node2int_file, 'w') as out:
            out.write(''.join(["%s\t%s\n" % (n, node2int[n]) for n in node2int]))

    #goid_prots = parse_gain_file(gain_file)
    # get the positives and negatives from the matrix
    #pos_neg_file = "inputs/pos-neg/all/pos-neg-bp-100.tsv"
    pos_neg_file = "inputs/pos-neg/noniea/nonieapos-neg-bp-100.tsv"
    goid_prots, goid_neg = parse_pos_neg_matrix(pos_neg_file)
    # TODO test the difference between the GAIN pos/neg and our method of assigning positives and negatives
    #goterms = set(["9405"])
    goterms = set(["GO:0009405"])  # pathogenesis

    k = 200
    t = 2
    #s = 20
    #s_list = [10, 25, 50]
    #s_list = [50, 200]
    #s_list = [20]
    s_list = [200]
    #s_list = [500, 100]
    #s_list = range(2,10,2)
    #s_list = [100, 200, 500, 750, 1000][::-1]
    a = 1
    #a = 0.8
    eps = 0.001
    #eps = None

    # make a copy of the graph to try using SinkSource's UB and LB on SinkSourcePlus
    # Add a fake negative (ground)?
    # Or I could just use alpha and normally set it =1 for SinkSource!

    for goterm in goterms:
        positives = set([node2int[n] for n in goid_prots[goterm] if n in node2int])
        negatives = set([node2int[n] for n in goid_neg[goterm] if n in node2int])
        positives = positives.intersection(set(H.nodes()))
        negatives = negatives.intersection(set(H.nodes()))
        # covert the positives to integers to match the graph

    #    for k in 
#        for alg in algorithms:
#            # This uses the same UB as Ripple
#            if alg == "sinksourceplus-ripple":
#            # This uses the same UB as Ripple, but with negatives
#            if alg == "sinksource-ripple":
#            # This computes UB and LB using iteration
#            if alg == "sinksource-topk-ub":
#            # This computes UB and LB using iteration, but doesn't use negatives
#            if alg == "sinksourceplus-topk-ub":
        topk_node_scores_by_s = {}
        topk_node_times_by_s = {}
        topk_node_iters_by_s = {}
        for s in s_list:
            print("Running Top-k SinkSourceRipple with %d positives for GO term %s" % (len(positives), goterm))
            #R, scores, time, iters = sinksource_topk_ub.runSinkSourceTopK(H, positives, negatives, k=k, t=t, s=s, eps=eps)
            #R, scores, time, iters = sinksource_topk_ub.runSinkSourceTopK(H, positives, k=k, t=t, s=s, a=a, eps=eps)
            # TODO also try adding a alpha to sinksource
            R, scores, time, iters = sinksource_topk_ub.runSinkSourceTopK(H, positives, negatives=negatives, k=k, t=t, s=s, a=a, eps=eps)
            #R, scores, time, iters = sinksource_ripple.runSinkSourceRipple(H, positives, k=k, t=t, s=s, a=a, eps=eps)
            # TODO also try adding negatives to sinksourceplus
            #R, scores, time, iters = sinksource_ripple.runSinkSourceRipple(H, positives, negatives=negatives, k=k, t=t, s=s, a=a, eps=eps)
            #node_names = int2nodescores
            # make sure the topk remain the same
            topk_node_scores_by_s[s] = {n: scores[node2int[n]] for n in sorted([int2node[u] for u in R])}
            topk_node_times_by_s[s] = time
            topk_node_iters_by_s[s] = iters
        for s in s_list:
            scores = topk_node_scores_by_s[s]
            print("s=%d: " % (s) + ', '.join(["%s: %0.3f" % (n, scores[n]) for n in sorted(scores)]))
        # print out a table of the times
        print("k\tt\ta\ts\tsec\titers")
        #print("Times (sec):")
        for s in sorted(topk_node_times_by_s):
            print("%d\t%d\t%s\t%d\t%0.2f\t%d" % (k, t, str(a), s, topk_node_times_by_s[s], topk_node_iters_by_s[s]))
        #print(', '.join(["%d: %0.2f" % (s, topk_node_times_by_s[s]) for s in sorted(topk_node_times_by_s)]))


    #    k = 50
    #    print("Running SinkSourcePlus with %d positives for GO term %s" % (len(positives), goterm))
    #    out_file = "test/sinksource/node-scores-ub-k%d-%s-%s-e0_0001.txt" % (k, version, goterm)
    #    if os.path.isfile(out_file):
    #        print("%s already exists. Skipping running regular SinkSourcePlus")
    #        lines = utils.readColumns(out_file, 1, 2)
    #        scores_ub = {n: s for n,s in lines}
    #        sort_scores_ub = [n for n,s in lines]
    #    else:
    #        sinksourceplusG_ub = sinksource.runSinkSource(G.copy(), positives, k=k)
    #        scores_ub = nx.get_node_attributes(sinksourceplusG_ub, 's')
    #        sort_scores_ub = sorted(scores_ub, key=scores_ub.get, reverse=True)
    #        print("Writing scores of top %d nodes to %s" % (k, out_file))
    #        # write it to a file
    #        with open(out_file, 'w') as out:
    #            out.write(''.join(["%s\t%0.3f\n" % (n, scores_ub[n]) for n in sort_scores_ub]))
    #
    #    out_file = "test/sinksource/node-scores-%s-%s-e0_00001.txt" % (version, goterm)
    #    if os.path.isfile(out_file):
    #        print("%s already exists. Skipping running regular SinkSourcePlus")
    #        lines = utils.readColumns(out_file, 1, 2)
    #        scores = {n: s for n,s in lines}
    #        sort_scores = [n for n,s in lines]
    #    else:
    #        sinksourceplusG = sinksource.runSinkSource(G.copy(), positives)
    #
    #        scores = nx.get_node_attributes(sinksourceplusG, 's')
    #        sort_scores = sorted(scores, key=scores.get, reverse=True)
    #        print("Writing scores to %s" % (out_file))
    #        # write it to a file
    #        with open(out_file, 'w') as out:
    #            out.write(''.join(["%s\t%0.3f\n" % (n, scores[n]) for n in sort_scores]))
    #
    #    #scores_ub = nx.get_node_attributes(sinksourceplusG_ub, 's')
    #    #sort_scores_ub = sorted(scores_ub, key=scores_ub.get, reverse=True)
    #    print("node\tUB-score\tReg-score\tReg-score-idx")
    #    for n in sort_scores_ub[:50]:
    #        reg_idx = sort_scores.index(n)
    #        print("%s\t%0.2f\t%0.2f\t%d" % (n, scores_ub[n], scores[n], reg_idx))
    #    print("node\tReg-score\tUB-score\tUB-score-idx")
    #    for n in sort_scores[:50]:
    #        ub_idx = sort_scores_ub.index(n)
    #        print("%s\t%0.2f\t%0.2f\t%d" % (n, scores[n], scores_ub[n], ub_idx))


def parse_args(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--version',type='string', action='append',
                      help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    parser.add_option('', '--algorithm', type='string', metavar='STR',
                      help="Algorithm for which to get predictions. Options: '%s'" % ("', '".join(f_settings.ALGORITHM_OPTIONS.keys())))
    parser.add_option('', '--exp-name', type='string',
                      help="Experiment name to use when running GAIN.")
    parser.add_option('', '--goid', type='string', metavar='STR',
                      help='GO-term ID for which annotations and precitions will be posted')

    (opts, args) = parser.parse_args(args)

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)
    if opts.exp_name is None or opts.goid is None or opts.algorithm is None:
        print("--exp-name, --goid, --algorithm, required")
        sys.exit(1)

    #if opts.algorithm not in f_settings.ALGORITHM_OPTIONS:
    #    print "--algorithm %s is not a valid algorithm name" % (opts.algorithm)
    #    sys.exit(1)

    # TODO
    for alg in opts.algorithm:
        if alg not in ALGORITHMS:
            print("ERROR: '%s' not a valid algorithm name. Algorithms are: '%s'." % (alg, ', '.join(ALGORITHMS)))
            sys.exit(1)

    return opts


if __name__ == "__main__":
    main()
