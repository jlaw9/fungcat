#!/usr/bin/python

print("Importing libraries")

import os
import sys
import glob
import shutil
import itertools
from os import path
from optparse import OptionParser
# f_settings is in the folder above. Use this to import it
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
import utils.file_utils as utils
import fungcat_settings as f_settings
import matplotlib 
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
# for the shape version 
import matplotlib.lines as mlines
import pandas as pd
import seaborn as sns
sns.set_style('whitegrid')
from collections import defaultdict
from sklearn import metrics
from scipy.stats import kruskal, mannwhitneyu
# Use the functions and such from here:
import plot_prec_rec


def main(versions, exp_names, algorithms, selected_terms, test_sig=False, pdf=False, forced=False):
     
    version_cv_results, goids_info = plot_prec_rec.get_cv_results(versions, exp_names, algorithms, selected_terms)

    # compare the fmax, auprc, and auroc of the algorithms run
    # to see if there are any algorithms that significantly outperform the others.
    # First run the 

    # organize the figures by the exp names being combined
    out_dir = "outputs/viz/eval/%s" % ('-'.join(exp_names))
    utils.checkDir(out_dir)

    # make a boxplot of each of the figures


    alg_results = get_alg_measure_results(version_cv_results, goids_info, algorithms) 
    boxplots_all_results(
            alg_results, goids_info, algorithms, 
            out_pref='%s/' % (out_dir), pdf=pdf, forced=forced)

    #measure_name = {
    #        'fmax': "F-max",
    #        "auprc": "Area under Precision-Recall Curve",
    #        "avgp": "Weighted Average Precision",
    #        "auroc": "Area under ROC",
    #        }
    if opts.test_sig is True:
        # now use the Kruskal-Wallis test to test if the null hypothesis 
        # that the population median of all of the groups is equal
        for measure in ('fmax', 'auprc', 'avgp', 'auroc'):
        #for measure in [('fmax')]:
            version_alg_measure_results = [] 
            for version in sorted(versions):
                for alg in algorithms:
                    version_alg_measure_results.append(alg_results[version][alg][measure])
            print("Running kruskal-wallis test on %s results comparing %d boxplots" % (measure, len(version_alg_measure_results)))
            print(kruskal(*version_alg_measure_results))

        for measure in ('fmax', 'auprc', 'avgp', 'auroc'):
            # now compare the individual boxplots with the mann whitney rank test
            #for (v1, a1), (v2, a2) in itertools.combinations(list(itertools.product(versions, algorithms)), 2):
            print("Running Mann-Whitney U test for '%s' on each version combination and algorithm combination" % (measure))
            for a in algorithms:
                for v1, v2 in itertools.combinations(versions, 2):
                    #print("\tMann-Whitney U test comparing %s across versions %s with %s" % (alg, v1, v2))
                    #print(mannwhitneyu(alg_results[v1][alg][measure], alg_results[v2][alg][measure]))
                    test_statistic, pval = mannwhitneyu(alg_results[v1][a][measure], alg_results[v2][a][measure]) 
                    print("%s\t%s\t%s\t%s\t%0.3e" % (v1, v2, a, a, pval))
            for v in versions:
                for a1, a2 in itertools.combinations(algorithms, 2):
                    #print("\tMann-Whitney U test comparing algorithms %s - %s for %s" % (a1, a2, v))
                    #print(mannwhitneyu(alg_results[v][a1][measure], alg_results[v][a2][measure]))
                    test_statistic, pval = mannwhitneyu(alg_results[v][a1][measure], alg_results[v][a2][measure]) 
                    print("%s\t%s\t%s\t%s\t%0.3e" % (v, v, a1, a2, pval))



def compute_fmax(alg_prec_rec):
    f_measures = []
    for p,r in alg_prec_rec:
        harmonic_mean = (2*p*r)/(p+r)
        f_measures.append(harmonic_mean)
    return max(f_measures)


def compute_avgp(alg_prec_rec):
    # average precision score
    # see http://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html#sklearn.metrics.average_precision_score
    avgp = 0
    prev_r = 0 
    for p,r in alg_prec_rec:
        recall_change = r - prev_r
        avgp += (recall_change*p)
        prev_r = r
    #avgp = avgp / float(len(alg_prec_rec))
    return avgp


def compute_auprc(alg_prec_rec):
    x = [r for p,r in alg_prec_rec]
    y = [p for p,r in alg_prec_rec]
    auprc = metrics.auc(x, y, reorder=True)
    return auprc


def compute_auroc(alg_tpr_fpr):
    x = [fpr for tpr,fpr in alg_tpr_fpr]
    y = [tpr for tpr,fpr in alg_tpr_fpr]
    auroc = metrics.auc(x, y)
    return auroc


#def plot_boxplots():
def boxplots_all_results(
        alg_results, goids_info,
        algorithms, title="Scatterplot of F-max results for SinkSource",
        out_pref=None, pdf=False, forced=False):

    measure_name = {
            'fmax': "F-max",
            "auprc": "Area under Precision-Recall Curve",
            "avgp": "Weighted Average Precision",
            "auroc": "Area under ROC",
            }
    for measure in ('fmax', 'auprc', 'avgp', 'auroc'):
        print("Plotting boxplots for %s" % (measure))
        out_file = "%sboxplots-%s-%s.png" % (out_pref, measure, '-'.join(sorted(alg_results.keys())))

        # number of rows is the number of versions
        fig, axes = plt.subplots(nrows=len(alg_results), figsize=(6,8))
        if len(alg_results) == 1:
            axes = [axes]
        plt.suptitle("Boxplots of %s results\nfor %d GO terms" % (measure_name[measure], len(goids_info)))
        for version, ax in zip(sorted(alg_results.keys()), axes):
            alg_measure_results = [alg_results[version][alg][measure] for alg in algorithms]
            bplot = ax.boxplot(alg_measure_results, 
                    labels=[plot_prec_rec.NAMES[alg] for alg in algorithms])
                    #patch_artist=True)  # fill with color
            #sns.boxplot(y=alg_measure_results,
            #        labels=[plot_prec_rec.NAMES[alg] for alg in algorithms], 
            #        ax=ax)

            # taken from here: https://matplotlib.org/examples/statistics/boxplot_color_demo.html
            #for patch, color in zip(bplot['boxes'], [plot_prec_rec.COLORS[alg] for alg in algorithms]):
            #    patch.set_facecolor(color)
            ax.set_title(version)
            ax.set_xticklabels([plot_prec_rec.NAMES[alg] for alg in algorithms], rotation=10)
            ax.set_ylim(ymax=1.00, ymin=-0.01)
        #legend = ax.legend(frameon=True)
        ## Put a nicer background color on the legend.
        #legend.get_frame().set_facecolor('#ffffff')
        #legend.get_frame().set_linewidth(2.0)
        ## if the title is too long, the plot does weird things
        #plt.tight_layout()

        plot_prec_rec.savefig(fig, out_file, pdf=pdf, forced=forced)


def get_alg_measure_results(version_cv_results, goids_info, algorithms):
    # get the fmax of the given algorithm for each version
    alg_results = {version: {alg: defaultdict(list) for alg in algorithms} for version in version_cv_results}
    for version in version_cv_results:
        for goid in version_cv_results[version]:
            for alg in version_cv_results[version][goid]:
                alg_cv = version_cv_results[version][goid][alg]
                if len(alg_cv) == 0 or len(alg_cv) == 1:
                    print(version, goid, alg, "has <= 1 CV line. Num genes: %d. Skipping" % (int(goids_info[goid]['genes'])))
                    continue
                alg_prec_rec = [(p,r) for p,r,fpr in alg_cv]
                # keep the recall (TPR) and false positive rate (FPR)
                alg_tpr_fpr = [(r,fpr) for p,r,fpr in alg_cv]
                alg_results[version][alg]['fmax'].append(compute_fmax(alg_prec_rec))
                alg_results[version][alg]['auprc'].append(compute_auprc(alg_prec_rec))
                alg_results[version][alg]['avgp'].append(compute_avgp(alg_prec_rec))
                alg_results[version][alg]['auroc'].append(compute_auroc(alg_tpr_fpr))

    return alg_results


#def plot_boxplots():
if __name__ == '__main__':
    opts, args = plot_prec_rec.parseArgs(sys.argv)

    algorithms = opts.algorithm
    if algorithms is None:
        algorithms = [
            "local-ova",
            "local-ovn",
            "fun-flow-ovn",
            "sinksource-ova",
            "sinksource-ovn",
            "genemania-ova",
        ]
    #for version in opts.version:
    only_functions = set() 
    if opts.only_functions:
        for f in opts.only_functions:
            functions = utils.readItemSet(f, 1)
            print("read %d functions from %s" % (len(functions), f))
            only_functions.update(functions)
    else:
        only_functions.add(opts.selected_term) 
    main(opts.version, opts.exp_name, algorithms, only_functions, test_sig=opts.test_sig, pdf=opts.pdf, forced=opts.forced)
