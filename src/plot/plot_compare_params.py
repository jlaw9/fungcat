# Compare parameters for each algorithm

from optparse import OptionParser, OptionGroup
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
cutoffs = [(10,1000)]


def main(versions, **kwargs):
    #kwargs['main_weight_type'] = plt_utils.MAIN_WEIGHT_TYPE 
    #kwargs['main_alg_opts'] = plt_utils.MAIN_ALG_OPTS 
    #kwargs['alg_names'] = plt_utils.ALG_NAMES 
    if len(versions) > 1 or len(kwargs['algorithm']) > 1 or \
       len(kwargs['ev_codes']) > 1 or len(kwargs['eval_ev_codes']) > 1:
        print("\nError: len(versions) > 1 or len(algorithms) > 1")
        print("This script is for comparing the params of 1 version, 1 alg at a time")
        print("Quitting.")
        sys.exit()

    ev_code_results = plt_utils.get_results(versions, **kwargs)
    measures = ["fmax"]
    plt_utils.results_overview(ev_code_results, measures=measures, **kwargs) 

    alg = kwargs['algorithm'][0] 
    # pull out the single version
    version, exp_name = list(ev_code_results.keys())[0]
    df = ev_code_results[(version, exp_name)]

    # make sure the # annotations is above 10
    # TODO make an option for this
    df = df[df['# test ann'] >= 10] 

    # pull out the fmax and reformat the dataframe
    print(sorted(df['Algorithm'].unique()))
    #for alg_params in df_curr['Algorithm'].unique():

    # reshape the dataframe to get a column per algorithm
    df_fmax = df.set_index(['#taxon', 'goid'])
    df_fmax = df_fmax[['fmax', 'Algorithm']]
    df_fmax = df_fmax.pivot(columns='Algorithm')
    #print(df.head(2))
    #pdb.set_trace()

    ## read the files
    #params_fmax = parse_files(
    #    versions, algorithms, exp_name, l_list, 
    #    ['-a'+str(a).replace('.','_') for a in a_list], eps_list, 
    ##     maxi_list, unweighted,
    #    ['-maxi'+str(m).replace('.','_') for m in maxi_list], unweighted,
    #    cv_file_template=cv_file_template, loso=True)
    colors = "Set2"
    params_list = []
    if len(kwargs['max_iters']) > 1 and len(kwargs['alpha']) > 1:
        ylabel = "Max # Iterations & Alpha"
        # if comparing both alpha and maxi, use this:
        params_list = ['a%s-%d' % (col.split(' a')[-1].split('-maxi')[0], int(col.split('-maxi')[-1])) for fmax, col in df_fmax.columns]
        exp_type = "maxi-alpha"
    elif len(kwargs['max_iters']) > 1:
        colors = "Set3"
        ylabel = "Max # Iterations"
        exp_type = "maxi-a%s" % (str(kwargs['alpha'][0]).replace('.','_'))
        #params_list = [int(col.split('-maxi')[-1]) for col in df_fmax.columns]
        params_list = kwargs['max_iters']
    elif len(kwargs['alpha']) > 1 and alg == 'sinksource':
        ylabel = "Alpha"
        exp_type = "alpha"
        params_list = kwargs['alpha']
    #elif algorithms[0] == 'genemania':
    elif len(kwargs['tol']) > 1:
        ylabel = "Tol"
        exp_type = "tol"
        #params_list = [float(col.split('-maxi20-tol')[-1].replace('_','.')) for col in df_fmax.columns]
        params_list = kwargs['tol']
        # TODO
    # also comapre birgrank params
    elif alg in ['birgrank', 'aptrank']:
        exp_type = ''
        for param, param_str in [('alpha', 'a'), ('eps', 'eps'), ('theta', 't'), ('mu', 'm'), ('br_lambda', 'l')]:
            if len(kwargs[param]) > 1:
                ylabel = param
                exp_type = "%s%s" % (param, exp_type)
                params_list = kwargs[param]
            else:
                exp_type += "-%s%s" % (param_str, str(kwargs[param][0]).replace('.','_'))

    #     a_list = [int(col.split('-maxi')[-1]) for col in df_fmax.columns]
    # make sure the columns match up with the dataframe, which sorts by string
    cols = sorted(str(x) for x in params_list)
    print("Exp_type: %s; renaming columns to: %s" % (exp_type, cols))
    if exp_type == "maxi-alpha":
        df_fmax.columns = cols
    else:
        df_fmax.columns = [float(x) for x in cols]
    print(df_fmax.head(2))

    ev_codes, eval_ev_codes, weight_type = kwargs['ev_codes'][0], kwargs['eval_ev_codes'][0], kwargs['weight_type'][0]
    eval_str = "-recov-"+eval_ev_codes if len(eval_ev_codes) > 1 else eval_ev_codes
    weight_str = '-'+weight_type if 'string' in version else ''
    out_dir = "outputs/viz/eval-loso/%s%s/%s/%s" % (ev_codes, eval_str, version, alg)
    utils.checkDir(out_dir)
    out_file = "%s/%s%s-compare-%s%s.pdf" % (out_dir, exp_name, weight_str, exp_type, kwargs['postfix'])

    title = "LOSO comparing %s %s \n %s, %s \n %d GO terms with >= 10 ann" % (
        exp_type, alg, version, exp_name, 
        len(set(df_fmax.index.get_level_values(1))))
    #if kwargs['line_plot']:
    if True:
        plot_loso_comparison_line(df_fmax, params_list, out_file=out_file, title=title, ylabel=ylabel, color_palette=colors, forced=kwargs['forceplot'])
    else:
        plot_loso_comparison(df_fmax, params_list, out_file=out_file, title=title, ylabel=ylabel, color_palette=colors, forced=kwargs['forceplot'])

    #print("medians: %s" % (', '.join(["%s: %0.3f" % (str(param), df_fmax[param].median()) for param in params_list])))
    # also write the medians to a file
    stat_file = out_file.replace('.pdf', '.txt')

    compute_stats(df_fmax, params_list, stat_file, alg=alg, exp_type=exp_type, 
            forced=kwargs['forceplot'], compare_to=kwargs['compare_to'])


def plot_loso_comparison(df_fmax, params_list, out_file='', title='', ylabel="Alpha", color_palette="Set3", forced=False):
    if forced or not os.path.isfile(out_file):
        print("Writing to %s" % (out_file))
    else:
        print("File already exists %s. Use --forceplot to overwrite it" % (out_file))
        return
    fig, ax = plt.subplots(figsize=(5,4))
    sns.boxplot(data=df_fmax, orient="h", order=params_list, palette=color_palette, fliersize=2, ax=ax)

    ax.set_xlim(-0.02, 1.05)
    # sns.boxplot(data=df, orient="h")
    # plt.xlabel(r'$F_{\mathrm{max}}$', fontsize=12, fontweight='bold')
    plt.xlabel(r'$F_{\mathrm{max}}$', fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.title(title)
    # plt.title("LOSO comparing alpha (1000 max iters) \n %s, %s \n %d %s GO terms with >= 10 ann, %s" % (
    # plt.tight_layout()
    plt.savefig(out_file, bbox_inches="tight")
    plt.show()
    plt.close()


def plot_loso_comparison_line(df_fmax, params_list, out_file='', title='', ylabel="Alpha", color_palette="Set3", forced=False):
    if forced or not os.path.isfile(out_file):
        print("Writing to %s" % (out_file))
    else:
        print("File already exists %s. Use --forceplot to overwrite it" % (out_file))
        return
    # melt the df back to full values
    df = df_fmax.melt(var_name='cols', value_name='vals')
    fig, ax = plt.subplots()
    #sns.boxplot(data=df_fmax, orient="h", order=params_list, palette=color_palette, fliersize=2, ax=ax)
    #sns.pointplot(x='cutoff', y='fmax', data=df_fmax, ax=ax
    sns.pointplot(x='cols', y='vals', data=df, ax=ax,
#                  order=params_list, 
#                   hue='Algorithm', hue_order=[alg_name[a] for a in algorithms],
#                  palette=my_palette,
                  estimator=np.median, ci=None, 
                #markers="o", color="#5f9e6e",
                markers="^", color="#b55d60",
                #markers="x", color="#857aaa",
                #markers="s", color="#c1b37f",
                #markers="P", color="#71aec0",
                 )

    plt.setp(ax.lines,linewidth=1)  # set lw for all lines of g axes

    ax.set_ylim(0.35, 0.475)
    ax.set_ylabel(r'$F_{\mathrm{max}}$', fontsize=12, fontweight='bold')
    ax.set_xlabel(ylabel)
    #ax.set_ylim(0.3, 0.625)
    plt.title(title)
#     plt.xticks(rotation=45)
    plt.savefig(out_file)
    plt.show()
    plt.close()


def compute_stats(df_fmax, params_list, stat_file, alg='sinksource', exp_type='alpha', forced=False, compare_to=None):
    if forced or not os.path.isfile(stat_file):
        if compare_to is None:
            alg1 = max(params_list)
            if exp_type == "maxi-alpha":
                alg1 = "a1.0-1000" 
        else:
            alg1 = compare_to
        alpha_1 = df_fmax[alg1].dropna().values
        out_str = "#alg\t%s\t%s-2\t2-median fmax\tpval\tmedian of differences\n" % (exp_type, exp_type)
        # for a in params_list[:-1]:
        for a in params_list:
            fmax_a = df_fmax[a].dropna().values
            test_statistic, pval = mannwhitneyu(alpha_1, fmax_a, alternative='greater') 
            out_str += "%s\t%s\t%s\t%0.3f\t%0.3e\t%0.3e\n" % (alg, max(params_list), str(a), np.median(fmax_a), pval, 
                    np.median(alpha_1 - fmax_a) if len(alpha_1) == len(fmax_a) else -1)
        print("appending to %s" % (stat_file))
        with open(stat_file, 'w') as out:
            out.write(out_str)
        print(out_str)
    else:
        print("%s already exists. Skipping" % (stat_file))


def parse_args(args):
    parser = plt_utils.setup_optparse()
    plt_utils.add_plot_opts(parser)

    # also add a couple of options specific to this script
    group = OptionGroup(parser, 'Stat-sig options')
    group.add_option('', '--compare-to', 
            help="Specify which parameter value to compare against when computing p-values. For example: 0.95-10. Default is the max param value (or a1.0-1000)")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args(args)
    kwargs = plt_utils.validate_opts(opts) 

    return kwargs
    

if __name__ == "__main__":
    kwargs = parse_args(sys.argv)
    versions = kwargs.pop('version')
    main(versions, **kwargs)
