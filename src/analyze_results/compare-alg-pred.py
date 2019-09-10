#!/usr/bin/env python

import sys, os
from os import path
from optparse import OptionParser
# f_settings is in the folder above. Use this to import it
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import fungcat_settings as f_settings
import utils.file_utils as utils
import fungcat_utils as f_utils

#import matplotlib 
#matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
#import matplotlib.pyplot as plt
## for the shape version 
#import matplotlib.lines as mlines
#import pandas as pd
#import seaborn as sns
#sns.set_style('whitegrid')
from collections import defaultdict
import numpy as np
from scipy.stats import describe


# for each GO term, compare how similar the predictions of two algorithms (or parameters) are
def main(files_to_compare):
    print("Comparing %d files" % (len(files_to_compare)))

    parsed_files = {}
    for i, f1 in enumerate(files_to_compare):
        goterm_preds1 = parse_goterm_pred(f1, parsed_files)
        for j, f2 in enumerate(files_to_compare):
            if i <= j:
                continue
            print("comparing %s with %s" % (f1, f2))
            goterm_preds2 = parse_goterm_pred(f2, parsed_files)

            goterms_in_common = goterm_preds1.keys() & goterm_preds2.keys()
            print("\t%d goterms in common" % (len(goterms_in_common)))
            goterm_jaccards = {}
            for goterm in goterms_in_common:
                 p1 = goterm_preds1[goterm]
                 p2 = goterm_preds2[goterm]
                 jaccard = len(p1 & p2) / float(len(p1 | p2))
                 goterm_jaccards[goterm] = jaccard 
                 #print("jaccard for %s: %0.4f" % (goterm, jaccard))
            jaccards = list(goterm_jaccards.values())
            print("\tAvg jaccard overlap: %0.3f (+- %0.3f)" % (
                np.mean(jaccards), np.std(jaccards)))
            print("\tmin: %0.3f, max: %0.3f" % (min(jaccards), max(jaccards)))
            print(describe(jaccards))


def parse_goterm_pred(goterm_pred_file, parsed_files):
    #with open(goterm_pred_file, 'r') as f:
    #    for line in f
    if goterm_pred_file in parsed_files:
        goterm_preds =  parsed_files[goterm_pred_file]
    else:
        lines = utils.readColumns(goterm_pred_file, 1, 2)
        goterm_preds = defaultdict(set)
        for goterm, prot in lines:
            goterm_preds[goterm].add(prot) 

    return goterm_preds


def parse_args(args):
    ## Parse command line args.
    usage = '%s [options] <file1> <file2>\n' % (sys.argv[0]) + \
            '\n\tCompare predictions of two algorithms/parameters (output by my python gain)'  
    parser = OptionParser(usage=usage)

    (opts, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_help()
        sys.exit()

    return opts, args


if __name__ == '__main__':
    opts, args = parse_args(sys.argv)
    main(args)
