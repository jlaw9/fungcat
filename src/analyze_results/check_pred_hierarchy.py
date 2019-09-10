#! /usr/bin/python

# script to check if prediction scores for all GO terms are consistent with the hierarchy
# i.e., as you go up the DAG, the score of a given node should only increase or stay the same

# example call: 
# with swsn, the network is the same for each term, and the predictions are consistent
#python src/analyze_results/check_pred_hierarchy.py --version 2018_06-seq-sim-e0_1-string -W 0 --only-pred --string-core --weight-swsn -A sinksource --alpha 0.95 --eps 0.0 --max-iters 10  --exp-name test --pos-neg-file inputs/pos-neg/expc-rem-neg-comp-iea/pos-neg-bp-10-list.tsv -G GO:0009432 -G GO:0006974 -G GO:0033554 -T 208964
# for weight-per-goterm, most nodes are consistent (95%)
#python src/analyze_results/check_pred_hierarchy.py --version 2018_06-seq-sim-e0_1-string -W 0 --only-pred --string-core --weight-per-goterm -A sinksource --alpha 0.95 --eps 0.0 --max-iters 10  --exp-name test --pos-neg-file inputs/pos-neg/expc-rem-neg-comp-iea/pos-neg-bp-10-list.tsv -G GO:0009432 -G GO:0006974 -G GO:0033554 -T 208964


import sys, os
from collections import defaultdict
#from optparse import OptionParser
import networkx as nx
sys.path.append('src')
import utils.file_utils as utils
import fungcat_settings as f_settings
#import utils.graphspace.post_to_graphspace as gs
import algorithms.gain_scipy.alg_utils as alg_utils
import algorithms.gain_scipy.eval_leave_one_species_out as eval_sp
# to load the ontology graph
import igacat.go_term_prediction_examples.go_term_prediction_examples as pred_ex
#from pandas import read_csv
from scipy import sparse
import numpy as np
from tqdm import tqdm


def main(version, exp_name, taxon_goid_scores, goids, go_dag,
         taxon=None, unweighted="-goterm-weight",
         alpha="0_95", eps="0_0", maxi=10):

    if opts.non_pos_as_neg_eval is False:
        exp_name += "-use-neg" 
    else:
        exp_name += "-non-pos-neg" 

    # set the version for the UNIPROT_TO_SPECIES
    _, RESULTSPREFIX, NETWORK, _ = f_settings.set_version(version)

    ## for now, just use these. Make sure they are included in the options to the script
    #child       = "GO:0009432"  # SOS response
    #parent      = "GO:0006974"  # cellular response to DNA damage
    #grandparent = "GO:0033554"  # cellular response to stress
    goids2idx = {g: i for i, g in enumerate(goids)}
    #pred_results_file = "%s/%s/%s/loso%s-l0-a%s-eps%s-maxi%s.txt" % (
    #    RESULTSPREFIX, algorithm, exp_name, unweighted, alpha, eps, maxi)

    # now check the predictions to see if they are consistent with the hierarchy
    goid_perc_cons = {}
    #for t1, t2 in [(child, parent), (parent, grandparent), (child, grandparent)]:
    for t1 in goids:
        # only check terms with some non-negatie scores
        if taxon_goid_scores[goids2idx[t1]].nnz == 0:
            continue
        # get the parent terms for this term
        parent_terms = go_dag.neighbors(t1)
        for t2 in parent_terms:
            if t2 not in goids2idx or taxon_goid_scores[goids2idx[t2]].nnz == 0:
                continue
            t1_scores = taxon_goid_scores[goids2idx[t1]].toarray().flatten()
            t2_scores = taxon_goid_scores[goids2idx[t2]].toarray().flatten()
            # keep track of the % of nodes whose scores are consistent with the hierarchy
            num_consistent = 0
            for n, s1 in enumerate(t1_scores):
                #print(i, s1, t2_scores[i])
                if s1 <= t2_scores[n]:
                    num_consistent += 1
            perc_consistent = num_consistent / float(len(t1_scores))
            #print("%s, %s: %0.2f%% consistent" % (t1, t2, perc_consistent*100))
            goid_perc_cons[(t1, t2)] = perc_consistent

    return goid_perc_cons


if __name__ == '__main__':
    # I'm reusing the options from the eval_leave_one_species_out script
    opts, W, prots, ann_matrix, goids, test_ann_matrix = eval_sp.run(sys.argv)
    if 'bp' in opts.pos_neg_file:
        h = 'bp'
    else:
        h = 'mf'
    category = {'bp': 'P', 'mf': 'F'}
    # load the hierarchy to check parent-child relationships
    obo_file = f_settings.GO_FILE
    go_dags = pred_ex.parse_obo_file_and_build_dags(obo_file)
    print("Using %s GO DAG" % h.upper())
    go_dag = go_dags[category[h]]
    for goid in goids:
        if not go_dag.has_node(goid):
            print("ERROR: %s not in the %s GO DAG. Quitting." % (goid, h.upper()))
            sys.exit()

    # I need to re-run LOSO to get the prediction scores for each node and each GO term
    selected_species = utils.readDict(f_settings.VERSION_SELECTED_STRAINS[opts.version], 1, 2)
    taxons = opts.taxon
    if taxons is None:
        taxons = selected_species.keys()
    overall_perc_cons = {}  # dict of percent of nodes that are hierarchically consistent per goterm per taxon
    for taxon in tqdm(sorted(taxons)):
        print('\n'+'-'*50)
        print("Getting prediction scores for %d goterms for taxon %s" % (len(goids), str(opts.taxon)))
        taxon_goid_scores, params_results = eval_sp.main(opts.version, opts.exp_name, W, prots, ann_matrix, goids,
            opts.algorithm, opts, taxons=[taxon],
            eval_ann_matrix=test_ann_matrix)
        if taxon_goid_scores is None:
            continue

        goid_perc_cons = main(opts.version, opts.exp_name, taxon_goid_scores, goids, go_dag, taxon=taxon)
        print("\ttaxon: %s (%s); median %% nodes consistent for %d child-parent term pairs: %0.2f%%" % (
            selected_species[taxon], taxon, len(goid_perc_cons), np.median(list(goid_perc_cons.values()))*100))
        overall_perc_cons[taxon] = goid_perc_cons 

    all_perc = []
    print("species\ttaxon\t# term pairs\tmedian")
    for taxon in sorted(taxons):
        goid_perc_cons = overall_perc_cons[taxon]
        print("%s\t%s\t%d\t%0.2f%%" % (
            selected_species[taxon], taxon, len(goid_perc_cons), np.median(list(goid_perc_cons.values()))*100))
        all_perc += list(goid_perc_cons.values()) 
    print("all\t-\t%d\t%0.2f%%" % (len(all_perc), np.median(all_perc)*100))

    out_file = "outputs/stats/%s-%s-%s-%s-hier-consistent.txt" % (opts.version, opts.exp_name, h, opts.algorithm[0])
    print("writing to %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write(''.join("%s\t%s\t%s\t%0.2f\n" % (t, t1, t2, p*100) for t in sorted(overall_perc_cons) for (t1, t2), p in sorted(overall_perc_cons[t].items()))) 
