# Squeeze Evaluations - LOSO

from collections import defaultdict
import os
import sys 
from tqdm import tqdm
import itertools
sys.path.append("src")
import utils.file_utils as utils
import fungcat_settings as f_settings
import plot.plot_utils as plt_utils
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
# also compute the significance of sinksource vs local
from scipy.stats import kruskal, mannwhitneyu
import pdb


selected_species = ""
species_to_uniprot = ""
# TODO make an option for this
#cutoffs = [(10,1000)]


def main(**kwargs):
    #kwargs['main_weight_type'] = plt_utils.MAIN_WEIGHT_TYPE 
    #kwargs['main_alg_opts'] = plt_utils.MAIN_ALG_OPTS 
    #kwargs['alg_names'] = plt_utils.ALG_NAMES 
    if len(kwargs['version']) > 1 or len(kwargs['algorithm']) > 1 or \
       len(kwargs['ev_codes']) > 1 or len(kwargs['eval_ev_codes']) > 1:
        print("\nError: len(versions) > 1 or len(algorithms) > 1")
        print("This script is for plotting the SinkSourceSqueeze ranks comparison of a single version")
        print("Quitting.")
        sys.exit()

    alg = "sinksource-squeeze"
    version, ev_codes, eval_ev_codes, string_nets, h, pos_neg_str, weight_type, alpha = [
        kwargs[x][0] for x in ['version', 'ev_codes', 'eval_ev_codes', 'string_nets', 'hierarchy', 'pos_neg_str', 'weight_type', 'alpha']]
    size_range = kwargs['size_range']

    eval_str = "-recov-"+eval_ev_codes if len(eval_ev_codes) > 1 else eval_ev_codes
    weight_str = '-'+weight_type if 'string' in version else ''
    stats_dir = "outputs/%s/all/%s/%s-%s%s%s-%s-use-neg" % (
        version, alg, ev_codes, size_range, eval_ev_codes, '-'+string_nets if 'string' in version else '', h)
    out_dir = "outputs/viz/eval-loso/%s%s/%s/%s" % (ev_codes, eval_str, version, alg)
    #out_dir = "outputs/viz/ranks/loso-%s%s/%s" % (exp_name, eval_ev_codes, version)
    utils.checkDir(out_dir)

    # type of sinksource-squeeze rank comparison
    # exp_type = "all-"
    exp_type = "pos-neg-"
    # top-k for Squeeze
    k_str = 'all' 
    # this file is output by the --compare-ranks option for SS-Squeeze
    all_ranks_file = "%s/pred-compare-%sranks-a%s-k%s.txt" % (stats_dir, exp_type, str(alpha).replace('.','_'), k_str)
    # read files
    print("Reading rankings from %s" % (all_ranks_file))
    if not os.path.isfile(all_ranks_file):
        print("\tERROR: does not exist. Quitting")
        sys.exit()
    df = pd.read_csv(all_ranks_file, sep='\t')
    # keep only the GO terms with at least 10 annotations
    df = df[df['num_pos'] >= 10]
    # get all goterm-taxon pairs
    df['goterm-taxon'] = df['#goterm'] + '-' + df['taxon'].map(str)

    out_pref = "%s/%s%s-%s" % (out_dir, exp_type, h, alpha)
    print("out_pref: %s" % (out_pref))
    print("%d GO terms, %d taxon, %d GO term-taxon pairs" % (df['#goterm'].nunique(), df['taxon'].nunique(), df['goterm-taxon'].nunique()))
    print(df.head())
    #print(df[['#goterm', 'iter', 'kendalltau']].head(200))

    # for each goterm, get the iteration at which kendalltau hits 95%, 99% and 100%
    iter_80 = get_iteration_at_cutoff(df, 0.80, col_to_get='iter', cutoff_col='kendalltau', less_than=False)
    iter_90 = get_iteration_at_cutoff(df, 0.90, col_to_get='iter', cutoff_col='kendalltau', less_than=False)
    iter_95 = get_iteration_at_cutoff(df, 0.95, col_to_get='iter', cutoff_col='kendalltau', less_than=False)
    iter_99 = get_iteration_at_cutoff(df, 0.99, col_to_get='iter', cutoff_col='kendalltau', less_than=False)
    iter_100 = get_iteration_at_cutoff(df, 1.0, col_to_get='iter', cutoff_col='kendalltau', less_than=False)
    df_cutoffs = pd.DataFrame({'0.80': iter_80, '0.90': iter_90, '0.95': iter_95, '0.99': iter_99, '1.0': iter_100})
    df_cutoffs = df_cutoffs[['0.80', '0.90', '0.95', '0.99', '1.0']]
    #     - Also plot the total # of iterations it takes to fix all node ranks
    total_iters = df.groupby('goterm-taxon')['iter'].max()
    #total_iters.head()
    df_cutoffs['Fixed ordering'] = total_iters
    # for each goterm, get the iteration at which kendalltau hits 95%, 99% and 100%
    #df_cutoffs.head()
    for col in df_cutoffs.columns:
        print("%s median: %d (%d values)" % (col, df_cutoffs[col].median(), df_cutoffs[col].dropna().count()))
    #df_cutoffs.head()

    plot(df_cutoffs, out_pref)


def plot(df_cutoffs, out_pref):
    # fig, (ax1, ax2) = plt.subplots(ncols=2)

    # insert a break into the plot to better show the small and large ranges
    f, (ax2, ax1) = plt.subplots(ncols=1, nrows=2, sharex=True, figsize=(5,4.5))
    sns.boxplot(df_cutoffs, ax=ax1, order=["0.80", '0.90', '0.95', '0.99', '1.0', 'Fixed ordering'], fliersize=1.5)
    sns.boxplot(df_cutoffs, ax=ax2, order=["0.80", '0.90', '0.95', '0.99', '1.0', 'Fixed ordering'], fliersize=1.5)
    # f, ax1 = plt.subplots()
    # sns.boxplot(df_cutoffs, ax=ax1)

    ymin, ymax = ax1.get_ylim()
    # ax1.set_ylim(-2, df_cutoffs['1.0'].max())
    ax1.set_ylim(-2, df_cutoffs['0.99'].max()+5)
    #ax1.set_ylim(0, 75)
    ymin, ymax = ax2.get_ylim()
    ax2.set_ylim(df_cutoffs['0.99'].max()+10, ymax)
    #ax2.set_ylim(df_cutoffs['Fixed ordering'].min()-5, ymax)
    #ax2.set_ylim(df_cutoffs['1.0'].min()-10, ymax)

    # now get the fancy diagonal lines
    d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass to plot, just so we don't keep repeating them
    kwargs = dict(transform=ax2.transAxes, color='tab:gray', clip_on=False)
    # ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax2.plot((+d, -d), (-d, +d), **kwargs)        # top-left diagonal
    ax2.plot((1 + d, 1 - d), (-d, +d), **kwargs)  # top-right diagonal

    kwargs.update(transform=ax1.transAxes)  # switch to the bottom axes
    ax1.plot((+d, -d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 + d, 1 - d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    xlabel = r"Kendall's $\tau$ cutoff compared to fixed ordering"
    ylabel = "# Iterations"
    # ax1.set_xlabel(xlabel)
    # ax1.set_ylabel(ylabel)
    # ax1.set_title("# iterations to \n rank nodes correctly")
    # plt.suptitle("%d goterms, %s, %s k\n%s" % (df['#goterm'].nunique(), alg, k_str, exp_name))
    f.text(0.02, 0.5, ylabel, va='center', rotation='vertical', fontsize=12, fontweight='bold')
    f.text(0.5, 0.02, xlabel, ha='center', fontsize=12, fontweight='bold')
    # plt.tight_layout()

    # out_file = all_ranks_file.replace('.txt', '-%s.png' % (h))
    # out_file = "%s/%s-%s-loso-pos-neg-ktau-cutoffs-boxplots.png" % (out_dir, version, h)
    out_file = "%s-ktau-cutoffs-boxplots.pdf" % (out_pref)
    print("Writing figure to %s" % (out_file))
    plt.savefig(out_file)
    plt.show()
    plt.close()


def plot_series(s, title='', xlabel='', ylabel='', out_file=None):
    fig, ax = plt.subplots()
    s.index += 1
    s.plot()
    # also add an inlet 
    s.index -= 1
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if s.min() < 1e-10:
        ax.set_yscale('log')
    plt.tight_layout()
    plt.show()
    if out_file is not None:
        print("writing figure to %s" % (out_file))
        plt.savefig(out_file)
    plt.close()


def get_iteration_at_cutoff(df, cutoff, col_to_get='kendalltau', cutoff_col='max_d', less_than=True):
    """
    *less_than*: If True, find the first occurance <= cutoff. Otherwise find >= cutoff
    """
    val_at_cutoff = {}
    # UPDATE: should include both the goterm and taxon
    #for goterm_taxon in tqdm(df['goterm-taxon'].unique()):
    #    goterm_df = df[df['goterm-taxon'] == goterm_taxon]
    #    for v1, v2 in goterm_df[[col_to_get, cutoff_col]].values:
    goterm_taxons = df['goterm-taxon'].unique()
    df = df[['goterm-taxon', col_to_get, cutoff_col]]
    g = df.groupby(['goterm-taxon'])
    for goterm_taxon in tqdm(goterm_taxons):
        goterm_df = g.get_group(goterm_taxon)
        for v1, v2 in goterm_df[[col_to_get, cutoff_col]].values:
            if (less_than is True and v2 <= cutoff) \
               or (less_than is False and v2 >= cutoff):
                val_at_cutoff[goterm_taxon] = v1
                #print("%s: %s: %s; %s: %s" % (goterm, col_to_get, v1, cutoff_col, v2))
                #pdb.set_trace()
                #return
                break
        if goterm_taxon not in val_at_cutoff:
            val_at_cutoff[goterm_taxon] = v1
    #print(len(val_at_cutoff))
    #sys.exit()
    return val_at_cutoff

if __name__ == "__main__":
    kwargs = plt_utils.parse_args(sys.argv)
    main(**kwargs)
