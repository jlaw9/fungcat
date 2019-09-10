#! /usr/bin/python

# script to analyze the results of the s200 leave-one-species-out evaluation
# to see if the fact that IEA examples are removed as negatives
# is affecting the results for 200 species 


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
import algorithms.setup_sparse_networks as setup
#from pandas import read_csv
from scipy import sparse
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


def main(version, exp_name, prots, goterms=None):
    """
    """

#    if opts.non_pos_as_neg_eval is False:
#        exp_name += "-use-neg" 
#    else:
#        exp_name += "-non-pos-neg" 
#
#    # set the version for the UNIPROT_TO_SPECIES
#    _, RESULTSPREFIX, NETWORK, _ = f_settings.set_version(version)

    # TODO add options for all of these files and parameters
    # get the # of IEA ann for each GO term
    remnegiea_file = "inputs/pos-neg/2018_09/expc-comp-rem-neg-iea/pos-neg-bp-50-list.tsv"
    noniea_file = "inputs/pos-neg/2018_09/expc-comp/pos-neg-bp-50-list.tsv"
    #iea_file = "inputs/pos-neg/2018_09/iea/pos-neg-bp-50-list.tsv"
    ann_matrix_remnegiea, goids = setup.setup_sparse_annotations(remnegiea_file, goterms, prots,
                            selected_species=None, taxon=None)
    remnegiea_goids2idx = {g: i for i, g in enumerate(goids)}
    ann_matrix_noniea, goids = setup.setup_sparse_annotations(noniea_file, goterms, prots,
                            selected_species=None, taxon=None)
    noniea_goids2idx = {g: i for i, g in enumerate(goids)}
    # I can get the count from the fmax file and the summary stats file
    #ann_matrix_iea, iea_goids = setup.setup_sparse_annotations(iea_file, goterms, prots,
    #                        selected_species=None, taxon=None)
    #iea_goids2idx = {g: i for i, g in enumerate(goids)}
    iea_pos_neg_summary_file = "inputs/pos-neg/2018_09/iea/pos-neg-50-summary-stats.tsv"
    num_iea_ann = {goid: int(num_pos) for goid, num_pos in utils.readColumns(iea_pos_neg_summary_file, 1, 4) if goid != "GO term"}

    fmax_maxi1000_file = "outputs/2018_09-s200-seq-sim-e0_1/all/sinksource/expc-comp-rem-neg-iea-50-1000iea-bp-use-neg/loso-l0-a1_0-eps0_0-maxi1000.txt"
    fmax_maxi10_file = "outputs/2018_09-s200-seq-sim-e0_1/all/sinksource/expc-comp-rem-neg-iea-50-1000iea-bp-use-neg/loso-l0-a1_0-eps0_0-maxi10.txt"
    print("reading fmax for each sp-goterm pair from %s" % (fmax_maxi1000_file))
    fmax_maxi1000 = defaultdict(dict)
    all_goids = set()
    taxon_goid_num_iea_ann = {}
    with open(fmax_maxi1000_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            taxon, goid, fmax = line.rstrip().split('\t')[:3]
            all_goids.add(goid)
            taxon_goid_num_iea_ann[(taxon, goid)] = int(line.rstrip().split('\t')[-1])
            fmax_maxi1000[taxon][goid] = float(fmax) 
    print("reading fmax for each sp-goterm pair from %s" % (fmax_maxi10_file))
    fmax_maxi10 = defaultdict(dict)
    with open(fmax_maxi10_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            taxon, goid, fmax = line.rstrip().split('\t')[:3]
            fmax_maxi10[taxon][goid] = float(fmax) 

    print("Getting species of each prot from %s" % (f_settings.UNIPROT_TO_SPECIES))
    # for each of the 19 species, leave out their annotations 
    # and see how well we can retrieve them 
    uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1,2)
    # also build the reverse, but with the prot IDs
    node2idx = {n: i for i, n in enumerate(prots)}
    global species_to_uniprot_idx
    species_to_uniprot_idx = defaultdict(set)
    for p in uniprot_to_species:
        species_to_uniprot_idx[uniprot_to_species[p]].add(node2idx.get(p))
    for s in species_to_uniprot_idx:
        species_to_uniprot_idx[s].discard(None) 

    # now compute the difference in # negatives that have an IEA annotation per species
    # store a tuple of fmax diff and iea diff for each taxon - goid pair
    taxon_goid_ieadiff = {}
    fmax_diff = {}
    for goid in all_goids:
        if goid not in remnegiea_goids2idx:
            continue
        pos, neg = alg_utils.get_goid_pos_neg(ann_matrix_remnegiea, remnegiea_goids2idx[goid])
        remnegiea_neg = set(list(neg))
        pos, neg = alg_utils.get_goid_pos_neg(ann_matrix_noniea, noniea_goids2idx[goid])
        noniea_neg = set(list(neg))
        #pos, neg = alg_utils.get_goid_pos_neg(ann_matrix_iea, iea_goids2idx[goid])
        #iea_pos = set(list(pos))
        # TODO I should limit these to the proteins in the network(?)
        for taxon in fmax_maxi10:
            if goid not in fmax_maxi1000[taxon]:
                continue
            num_remnegiea_neg = len(remnegiea_neg - species_to_uniprot_idx[taxon])
            num_noniea_neg = len(noniea_neg - species_to_uniprot_idx[taxon])
            # need to normalize by something... How about the # of not-left-out IEA annotations
            #num_iea = len(iea_pos - species_to_uniprot_idx[taxon])
            num_iea = num_iea_ann[goid] - taxon_goid_num_iea_ann[(taxon, goid)]
            taxon_goid_ieadiff[(taxon, goid)] = float(num_noniea_neg - num_remnegiea_neg) / float(num_iea)
            fmax_diff[(taxon, goid)] = fmax_maxi1000[taxon][goid] - fmax_maxi10[taxon][goid]


    #print(taxon_goid_ieadiff)
    # now make a scatterplot comparing the iea diff with the fmax diff
    df = pd.DataFrame({'iea-diff': taxon_goid_ieadiff, 'fmax-diff': fmax_diff})
    bins = np.histogram(df['iea-diff'])[1]
    df['iea-group'] = pd.cut(df['iea-diff'], bins)
    df.boxplot('fmax-diff', by='iea-group')
    #grid = sns.jointplot(x="iea-diff", y="fmax-diff", data=df)
    plt.savefig("test2.png")
    plt.show()
    plt.close()


if __name__ == '__main__':
    # I'm reusing the options from the eval_leave_one_species_out script
    opts = eval_sp.parse_args(sys.argv)
    INPUTSPREFIX, _, net_file, selected_strains = f_settings.set_version(opts.version) 
    W, prots = alg_utils.setup_sparse_network(net_file)

    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 
    #goid_pos, goid_neg = alg_utils.parse_pos_neg_files([opts.pos_neg_file], goterms=goterms) 
    #_, goid_neg_noniea = alg_utils.parse_pos_neg_files([opts.pos_neg_file_eval], goterms=goterms) 
    # TODO

    main(opts.version, opts.exp_name, prots, goterms)
