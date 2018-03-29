
# Quick script to run/test the algorithms
print("Importing libraries")

from optparse import OptionParser
import os
import sys
from tqdm import tqdm
import itertools
import utils.file_utils as utils
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fungcat_settings as f_settings
import networkx as nx
import alg_utils
import sinksource
import sinksource_topk_ub
import sinksource_ripple
import sinksource_squeeze
from collections import defaultdict
#import pandas as pd


ALGORITHMS = [
    "sinksourceplus-ripple",  # This uses the same UB as Ripple
    "sinksource-ripple",  # This uses the same UB as Ripple, but with negatives
    "sinksource-topk-ub",  # This computes UB and LB using iteration
    "sinksourceplus-topk-ub",  # This uses the same UB as sinksource-ub-topk
    "sinksourceplus-squeeze",
    "sinksource-squeeze",
    "sinksourceplus",
    "sinksource",
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


def main(versions, exp_name, only_functions_file, pos_neg_file, algorithms=["sinksource-ripple"], eps=0.001, k_list=[100], t_list=[2], s_list=[50], a_list=[0.8], deltaUBLB_list=[None]):
    """
    *eps*: Convergence cutoff for sinksource and sinksourceplus
    """
    global int2node
    version_params_results = {}

    #version = "2017_10-seq-sim"
    #version = "2017_10-string"
    for version in versions:
        INPUTSPREFIX, RESULTSPREFIX, network_file, selected_strains = f_settings.set_version(version)
        #gain_file = f_settings.GOA_ALL_FUN_FILE
    #    # start with a single species
    #    #taxon = "208964"
    #    if taxon is not None:
    #        network_file = f_settings.STRING_TAXON_UNIPROT % (taxon, taxon, f_settings.STRING_CUTOFF)
    #        gain_file = f_settings.FUN_FILE % (taxon, taxon)

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
        #pos_neg_file = "inputs/pos-neg/noniea/nonieapos-neg-bp-100.tsv"
        goid_prots, goid_neg = parse_pos_neg_matrix(pos_neg_file)
        # organize the results by algorithm, then GO term, then algorithm parameters
        params_results = {}
        for alg in algorithms:
            out_dir = "%s/all/%s/%s" % (RESULTSPREFIX, alg, exp_name)
            utils.checkDir(out_dir)
            for goterm in sorted(goterms):
                if goterm not in goid_prots:
                    print("Skipping %s. Not in pos-neg-matrix file" % (goterm))
                    continue
                # covert the positives to integers to match the graph
                positives = set([node2int[n] for n in goid_prots[goterm] if n in node2int])
                negatives = set([node2int[n] for n in goid_neg[goterm] if n in node2int])
                positives = positives.intersection(set(H.nodes()))
                negatives = negatives.intersection(set(H.nodes()))
                out_pref = "%s/%s-" % (out_dir, goterm.replace(':',''))

                # now run the algorithm with all combinations of parameters
                curr_params_results = run_alg_with_params(H, alg, goterm, out_pref, positives, negatives,
                                                          k_list, t_list, s_list, a_list, deltaUBLB_list)
                params_results.update(curr_params_results)
                print(params_results_to_table(params_results))
        version_params_results[version] = params_results

    for version in versions:
        print(version)
        print(params_results_to_table(params_results))

    ## TODO test the difference between the GAIN pos/neg and our method of assigning positives and negatives
    ##goterms = set(["9405"])
    #goterms = set(["GO:0009405"])  # pathogenesis

#    k = 200
#    t = 2
#    #s = 20
#    #s_list = [10, 25, 50]
#    #s_list = [50, 100]
#    s_list = [200]
#    #s_list = [500, 1000, 2000]
#    #s_list = [500, 100]
#    #s_list = range(2,10,2)
#    #s_list = [100, 200, 500, 750, 1000][::-1]
#    #a = 1
#    a = 0.8
#    #deltaUBLB = 0.001
#    deltaUBLB = None

    # make a copy of the graph to try using SinkSource's UB and LB on SinkSourcePlus
    # Add a fake negative (ground)?
    # Or I could just use alpha and normally set it =1 for SinkSource!

def run_alg_with_params(H, alg, goterm, out_pref, positives, negatives,
                        k_list, t_list, s_list, a_list, deltaUBLB_list, forced=True):
    """
    """
    global int2node
    params_results = {}

    # TODO run these separately because they have different parameters
    if alg == "sinksourceplus" or alg == "sinksource":
        # TODO make eps a parameter
        eps_list = [0.0001, 0.001]
        for a, eps in tqdm(list(itertools.product(*[a_list, eps_list]))):
            out_file = "%sa%s-eps%s.txt" % (out_pref, str(a).replace('.', '_'), str(eps).replace('.', '_'))

            if forced is False and os.path.isfile(out_file):
                print("%s already exists. Skipping" % (out_file))
                continue

            if alg == "sinksourceplus":
                scores, time, iters = sinksource.runSinkSource(H, positives, negatives=None, max_iters=1000, delta=eps, a=a)
            if alg == "sinksource":
                scores, time, iters = sinksource.runSinkSource(H, positives, negatives=negatives, max_iters=1000, delta=eps, a=a)
            params_key = (alg, goterm, '-', '-', '-', a, 'eps=%s'%str(eps))
            params_results[params_key] = (time, iters, '-')

            print(params_results_to_table(params_results))
            print("Writing %s" % (out_file))
            # write the scores to a file
            # for now, write the top k node's results
            # TODO Possibly write the results to a JSON file so it's easier to load
            with open(out_file, 'w') as f:
                f.write("#%s\tGO term: %s\t%d positives\t%d negatives\ta=%s\teps=%s\n" \
                        % (alg, goterm, len(positives), len(negatives), str(a), str(eps)))
                for n in sorted(scores, key=scores.get, reverse=True):
                    f.write("%s\t%0.4f\n" % (int2node[n], scores[n]))
        return params_results

    # Run using all combinations of parameters
    # s and t won't change the results, just the amount of time, so loop through those separately
    total_params_list = list(itertools.product(*[k_list, t_list, s_list, a_list, deltaUBLB_list]))
    #params_list = list(itertools.product(*[k_list, a_list, deltaUBLB_list]))
    #t_s_list = list(itertools.product(*[t_list, s_list]))
    print("Running %d combinations of parameters" % (len(total_params_list)))
    for k, a, deltaUBLB in tqdm(list(itertools.product(*[k_list, a_list, deltaUBLB_list]))):

        #out_file = "%sk%d-a%s-d%s.txt" % (out_pref, k, str(a).replace('.', '_'),
        #                                    str(deltaUBLB).replace('.', '_'))
        #if os.path.isfile(out_file):
        #    print("%s already exists. Skipping running regular SinkSourcePlus" % (out_file))
        #else:
        for t, s in tqdm(list(itertools.product(*[t_list, s_list]))):
            print("Running %s with %d positives, %d negatives, k=%d, t=%d, s=%d, a=%s, deltaUBLB=%s for GO term %s" \
                    % (alg, len(positives), len(negatives), k, t, s, str(a), str(deltaUBLB), goterm))
            # TODO streamline calling the correct function. They all take the same parameters
            # This uses the same UB as Ripple
            if alg == "sinksourceplus-ripple":
                R, scores, time, iters, len_N = sinksource_ripple.runSinkSourceRipple(H, positives, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB)
            # This uses the same UB as Ripple, but with negatives
            elif alg == "sinksource-ripple":
                R, scores, time, iters, len_N = sinksource_ripple.runSinkSourceRipple(H, positives, negatives, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB)
            # This computes UB and LB using iteration
            elif alg == "sinksource-topk-ub":
                R, scores, time, iters, len_N = sinksource_topk_ub.runSinkSourceTopK(H, positives, negatives, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB)
            # This computes UB and LB using iteration, but doesn't use negatives
            elif alg == "sinksourceplus-topk-ub":
                R, scores, time, iters, len_N = sinksource_topk_ub.runSinkSourceTopK(H, positives, k=k, t=t, s=s, a=a, deltaUBLB=deltaUBLB)
            elif alg == "sinksourceplus-squeeze":
                R, scores, time, iters, len_N = sinksource_squeeze.runSinkSourceSqueeze(H, positives, k=k, a=a, deltaUBLB=deltaUBLB, ranked=True)
            elif alg == "sinksource-squeeze":
                R, scores, time, iters, len_N = sinksource_squeeze.runSinkSourceSqueeze(H, positives, negatives, k=k, a=a, deltaUBLB=deltaUBLB)


            out_file = "%sk%d-t%d-s%d-a%s-d%s.txt" % (out_pref, k, t, s, str(a).replace('.', '_'),
                                                str(deltaUBLB).replace('.', '_'))
            ## write the scores to a file
            ## for now, write the top X node's results
            #with open(out_file, 'w') as f:
            #    f.write("#%s\tGO term %s\t%d positives\t%d negatives\tk=%d\tt=%d\ts=%d\ta=%s\tdeltaUBLB=%s\n" % (goterm, alg, len(positives), len(negatives), k, s, str(a), str(deltaUBLB)))
            # Write the results to a JSON file so it's easier to load
            if os.path.isfile(out_file):
                print("%s already exists. Not re-writing to file" % (out_file))
            else:
                print("Writing %s" % (out_file))
                # write the scores to a file
                # for now, write the top k node's results
                # TODO Possibly write the results to a JSON file so it's easier to load
                with open(out_file, 'w') as f:
                    f.write("#%s\tGO term: %s\t%d positives\t%d negatives\tk=%d\tt=%d\ts=%d\ta=%s\tdeltaUBLB=%s\n" \
                            % (alg, goterm, len(positives), len(negatives), k, t, s, str(a), str(deltaUBLB)))
                    for n in sorted(R, key=scores.get, reverse=True):
                        f.write("%s\t%0.4f\n" % (int2node[n], scores[n]))

            # also keep track of the time it takes for each of the parameter sets
            params_key = (alg, goterm, k, t, s, a, deltaUBLB)
            params_results[params_key] = (time, iters, len_N)

            print(params_results_to_table(params_results))

    return params_results


def params_results_to_table(params_results):
    # print out a table of the results for each parameter set
    # print out the table after each iteration so the results are still viewable if the script quits earl or something
    results = "\t".join(['alg', 'goterm', 'k', 't', 's', 'a', 'deltaUBLB', 'time', 'iters', 'len_N']) + '\n'
    for params_key in sorted(params_results):
        alg, goterm, k, t, s, a, deltaUBLB = params_key
        time, iters, len_N = params_results[params_key]
        results += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.1f\t%d\t%s\n" % (alg, goterm, str(k), str(t), str(s), str(a), str(deltaUBLB), time, iters, str(len_N))
    return results

#        out_file = "test/sinksource/node-scores-%s-%s-%s-e0_001-2.txt" % (version, goterm.replace(':',''), str(a).replace('.', '_'))
#        if os.path.isfile(out_file):
#            print("%s already exists. Skipping running regular SinkSourcePlus" % (out_file))
#        else:
#            s = sinksource.runSinkSource(H, positives, negatives=None, max_iters=1000, delta=0.001, a=0.8)
#            sort_scores = sorted(s, key=s.get, reverse=True)
#            print("Writing scores to %s" % (out_file))
#            # write it to a file
#            with open(out_file, 'w') as out:
#                out.write(''.join(["%s\t%0.3f\n" % (int2node[n], s[n]) for n in sort_scores]))
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
    parser.add_option('-A', '--algorithm', action="append",
                      help="Algorithm for which to get predictions. Default is all of them. Options: '%s'" % ("', '".join(ALGORITHMS)))
    parser.add_option('', '--exp-name', type='string',
                      help="Experiment name to use when running GAIN.")
    parser.add_option('', '--pos-neg-file', type='string',
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
    parser.add_option('-d', '--deltaUBLB', type=float, action="append",
                      help="Parameter to fix node scores if UB - LB difference is < delta. Default=None")

    #parser.add_option('', '--goid', type='string', metavar='STR',
    #                  help='GO-term ID for which annotations and precitions will be posted')

    (opts, args) = parser.parse_args(args)

    if opts.exp_name is None or opts.pos_neg_file is None:
        print("--exp-name, --only_functions, --pos-neg-file, --algorithm, required")
        sys.exit(1)

    if opts.goterm is None and opts.only_functions is None:
        print("--goterm or --only_functions required")
        sys.exit(1)

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
    opts.deltaUBLB = opts.deltaUBLB if opts.deltaUBLB is not None else [None]

    return opts


if __name__ == "__main__":
    #versions = ["2017_10-seq-sim", "2017_10-seq-sim-x5-string"]
    opts = parse_args(sys.argv)

    goterms = set()
    if opts.only_functions is not None:
        only_functions = utils.readItemSet(opts.only_functions, 1)
        goterms = set(["GO:" + "0"*(7-len(str(x))) + str(x) for x in only_functions])
    if opts.goterm is not None:
        goterms.update(set(opts.goterm))

    main(opts.version, opts.exp_name, goterms,
            opts.pos_neg_file, opts.algorithm,
            k_list=opts.k, t_list=opts.t, s_list=opts.s, a_list=opts.alpha,
            deltaUBLB_list=opts.deltaUBLB)
    #main()
