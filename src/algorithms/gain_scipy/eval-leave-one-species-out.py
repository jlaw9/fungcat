
# script to see how well we can retrieve annotations of another species
# using the annotations of the other 18

from optparse import OptionParser,OptionGroup
from collections import defaultdict
import os
import sys
from tqdm import tqdm
import itertools
sys.path.append("src")
import utils.file_utils as utils
import fungcat_settings as f_settings
import run_algs

# start with 200 annotations, and then when processing the output, bin the fmax by the # of annotations per species
cut1 = 200
cut2 = 1000
#cut1 = 500
#cut2 = 1000
# use nothing for all evidence codes
evidence_codes = "all"
#evidence_codes = "rem-neg-iea"
#evidence_codes = "exp"
# TODO add command-line arguments
versions = ["2017_10-seq-sim"]
unweighted = False
#versions = ["2017_10-seq-sim-x5-string"]
#unweighted = True
# option to use negative examples when evaluating predictions
use_negatives_for_eval = False 
# start with less GO terms
#exp_name = "eval-species-interesting-terms" 
#only_functions_file = "inputs/only-functions/non-iea/interesting-terms.txt"
#exp_name = "species-rem-neg-iea-50"
exp_name = "eval-species-%s-%d-%d-25iters" % (evidence_codes, cut1, cut2)
if use_negatives_for_eval:
    exp_name += "-use-neg" 
only_functions_file = "inputs/only-functions/%s/%s-%d-%d.txt" % (evidence_codes, evidence_codes, cut1, cut2)
algorithms = ["sinksourceplus", "localplus", "sinksource"]
#algorithms = ["localplus", "sinksource"]
# selected species
# start with a single species and a single goterm
#selected_species = {"243277": ''}
selected_species_file = "inputs/selected-strains.txt"
selected_species = utils.readDict(selected_species_file, 1, 2)
#goterms = set(["GO:0009405"])
verbose = False

print("Getting species of each prot from %s" % (f_settings.UNIPROT_TO_SPECIES))
# for each of the 19 species, leave out their annotations 
# and see how well we can retrieve them 
uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1,2)
# also build the reverse
species_to_uniprot = defaultdict(set)
for p in uniprot_to_species:
    species_to_uniprot[uniprot_to_species[p]].add(p)

pos_neg_files = [
        "inputs/pos-neg/%s/pos-neg-bp-50-list.tsv" % (evidence_codes), 
        "inputs/pos-neg/%s/pos-neg-mf-50-list.tsv" % (evidence_codes)]
        #"inputs/pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-bp-50-list.tsv", 
        #"inputs/pos-neg/rem-neg-iea/rem-neg-iea-pos-neg-mf-50-list.tsv"]
goterms = run_algs.select_goterms(only_functions_file=only_functions_file) 
goid_pos, goid_neg = run_algs.parse_pos_neg_files(pos_neg_files, goterms=goterms) 
goid_pos = {goid: prots for goid, prots in goid_pos.items() if goid in goterms}
goid_neg = {goid: prots for goid, prots in goid_neg.items() if goid in goterms}

print("Training/testing with %d species, %d goterms" % (len(selected_species), len(goterms)))

for version in versions:
    for alg in algorithms:
        out_file = "outputs/%s/%s/%s/%s/ground-truth-l0-a0_8.txt" % (version, evidence_codes, alg, exp_name)
        if os.path.isfile(out_file):
            print("Removing %s as results will be appended to it for each taxon" % (out_file))
            os.remove(out_file)

# split the sets of positives/negatives into 19 sets with one species left out of each
for s in tqdm(sorted(selected_species)):
#for s in selected_species:
    # leave this taxon out by removing its annotations
    train_goid_pos = {}
    train_goid_neg = {}
    test_goid_pos = {}
    test_goid_neg = {}
    for goid in goid_pos:
        # TODO I should limit these to the proteins in the network
        train_prots = goid_pos[goid] - species_to_uniprot[s]
        test_prots = goid_pos[goid] & species_to_uniprot[s]
        train_neg = goid_neg[goid] - species_to_uniprot[s]
        test_neg = goid_neg[goid] & species_to_uniprot[s]
        if len(train_prots) == 0 or len(test_prots) == 0 or \
           (use_negatives_for_eval and (len(train_neg) == 0 or len(test_neg) == 0)):
            continue
        train_goid_pos[goid] = train_prots
        test_goid_pos[goid] = test_prots
        train_goid_neg[goid] = train_neg
        test_goid_neg[goid] = test_neg

    s_goterms = goterms & set(test_goid_pos.keys())
    tqdm.write("\n" + "-"*30)
    #print("\n" + "-"*30)
    #print("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
    tqdm.write("Taxon: %s - %s; %d/%d goterms with > 0 annotations" % (
        s, selected_species[s], len(s_goterms), len(goterms)))

    if len(s_goterms) == 0:
        print("\tskipping")
        continue

    if use_negatives_for_eval is False: 
        print("Evaluating using all non-ground-truth positives as false positives")
        test_goid_neg = None 
    else:
        print("Evaluating using only the ground-truth negatives labelled as positives as false positives")

    # for now, use most of the defaults
    # change to running sinksource with 25 iterations
    # leave alpha at default 0.8
    alg_runner = run_algs.Alg_Runner(
            versions, exp_name, train_goid_pos, train_goid_neg,
            s_goterms, algorithms=algorithms, eps_list=[0],
            unweighted=unweighted, max_iters=25,
            num_pred_to_write=0, verbose=verbose, 
            ground_truth=test_goid_pos, ground_truth_neg=test_goid_neg)
    # this will write an file containing the fmax for each goterm 
    # with the taxon name in the name of the file
    alg_runner.main(ground_truth_taxon=s)

