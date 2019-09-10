
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
import time


selected_species = ""
species_to_uniprot = ""
# TODO make an option for this
cutoffs = [(10,1000)]

def main(versions, **kwargs):
    kwargs['main_weight_type'] = plt_utils.MAIN_WEIGHT_TYPE 
    kwargs['main_alg_opts'] = plt_utils.MAIN_ALG_OPTS 
    kwargs['alg_names'] = plt_utils.ALG_NAMES 
    # alpha used for SS
    #kwargs['ss_alpha'] = plt_utils.MAIN_ALG_OPTS['sinksource'].split('-a')[-1].split('-eps')[0]

    global selected_species, species_to_uniprot
    selected_species = utils.readDict(f_settings.VERSION_SELECTED_STRAINS[versions[0]], 1, 2)
    uniprot_to_species = utils.readDict(f_settings.VERSION_UNIPROT_TO_SPECIES[versions[0]], 1,2)
    # also build the reverse
    species_to_uniprot = {s: set() for s in selected_species}
    for p in uniprot_to_species:
        species_to_uniprot[str(uniprot_to_species[p])].add(p)
    print(species_to_uniprot.keys())

    ev_code_results = get_results(versions, **kwargs)
    measures = ["fmax"]
    results_overview(ev_code_results, measures=measures, **kwargs) 
    generate_plots(ev_code_results, measures, **kwargs) 


def generate_plots(ev_code_results, measures=['fmax'], **kwargs):
    for version, ev_codes, eval_ev_codes, h in ev_code_results:
        cutoffs_data, out_dir = ev_code_results[(version, ev_codes, eval_ev_codes, h)]
        sort_taxon_by_fmax = None    

#         title = "Evaluation of recovery of %s ev. codes from %s ev. codes %s\n for 19 pathogenic bacteria. - %s  %s\n %d GO terms with %d <= # annotations < %d %s" % (
#             eval_ev_codes, ev_codes, use_neg_str, version, keep_ann, df_curr['goid'].nunique(), cutoff1, cutoff2, "overall" if split_overall else "")
        recov_str = " recovery of %s ev. codes from" % (eval_ev_codes) if eval_ev_codes != ""  else ""
        
        for i, (cutoff1, cutoff2) in enumerate(cutoffs):
            df_curr, species_exp_counts, species_comp_counts, species_iea_counts  = cutoffs_data[i]
            if 'comp' in eval_ev_codes:
                ann_stats = {'COMP': species_comp_counts}
            elif 'iea' in eval_ev_codes:
                ann_stats = {'ELEC': species_iea_counts}
            else:
                ann_stats = {'EXPC': species_exp_counts}
            title = ""
            if kwargs['for_pub'] is False:
                title = "Evaluation of%s %s ev. codes\n for 19 pathogenic bacteria. - %s\n %d %s GO terms with %d+ annotations" % (
                    recov_str, ev_codes, version, df_curr['goid'].nunique(), h.upper(), cutoff1)
            out_file = ""
            for measure in measures:
                print("Creating plots for '%s'" % (measure))

                weight_str = '-'+kwargs['weight_type'][0] if 'string' in version else ''
                out_file = "%s/%d-%d-%s-%s%s-%s%s%s.pdf" % (
                        out_dir, cutoff1, cutoff2, h, len(kwargs['algorithm']), kwargs['pos_neg_str'], measure, weight_str, kwargs['postfix'])
                plot_fmax_eval(df_curr, out_file, measure=measure, title=title,
                        sort_taxon_by_fmax=sort_taxon_by_fmax, ann_stats=ann_stats, **kwargs)

                out_file = out_file.replace('.pdf', '-sig.txt')
                # also get comparative p-values
                sig_results = eval_stat_sig(df_curr, out_file, measure=measure, sort_taxon_by_fmax=sort_taxon_by_fmax, **kwargs)

                # also write them to a file
                out_dir_stats = "%s/stats/" % (out_dir)
                utils.checkDir(out_dir_stats)
                out_file = "%s/%d-%d-%s-%s%s-%s%s-diff%s.tsv" % (
                    out_dir_stats, cutoff1, cutoff2, h, '-'.join([kwargs['alg1'].replace(' ','-'), kwargs['alg2'].replace(' ','-')]), kwargs['pos_neg_str'], measure, weight_str, kwargs['postfix'])
                scatterplot_fmax(df_curr, ev_codes, out_file, **kwargs)


def eval_stat_sig(
        df_curr, out_file, measure='fmax', 
        sort_taxon_by_fmax=None, alg1='sinksource', alg2='localplus', **kwargs):
    sp_pval = {}
    out_str = ""

    combinations = list(itertools.combinations(kwargs['algorithm'], 2))
    out_str += "#alg1\talg2\tpval\tCorrected p-value (x%d)\n" % (len(combinations))
    # Don't think sorting is needed
    #df_curr.sort_values('goid', inplace=True) 
    for a1, a2 in combinations:
        a1_fmax = df_curr[df_curr['Algorithm'] == kwargs['alg_names'].get(a1, a1)][measure]
        a2_fmax = df_curr[df_curr['Algorithm'] == kwargs['alg_names'].get(a2, a2)][measure]
        test_statistic, pval = mannwhitneyu(a1_fmax, a2_fmax, alternative='greater') 
        out_str += "%s\t%s\t%0.3e\t%0.3e\n" % (a1, a2, pval, pval*len(combinations))

    # also compare individual species
    curr_species = df_curr['#taxon'].unique()
    if sort_taxon_by_fmax is not None:
        curr_species = sort_taxon_by_fmax
    out_str += "Species\tAlg1\tAlg2\tAlg1-med\tAlg2-med\tRaw p-value\tCorrected p-value (x%d)\n" % (len(curr_species))
    for s in curr_species:
        name = f_settings.NAME_TO_SHORTNAME2.get(selected_species[str(s)],'-')
        df_s = df_curr[df_curr['#taxon'] == s]
        a1_fmax = df_s[df_s['Algorithm'] == kwargs['alg_names'].get(alg1, alg1)][measure]
        a2_fmax = df_s[df_s['Algorithm'] == kwargs['alg_names'].get(alg2, alg2)][measure]
        try:
            test_statistic, pval = mannwhitneyu(a1_fmax, a2_fmax, alternative='greater') 
            line = "%s\t%s\t%s\t%0.3f\t%0.3f\t%0.2e\t%0.2e" % (name, alg1, alg2, a1_fmax.median(), a2_fmax.median(), pval, pval*len(curr_species))
        except ValueError:
            line = "%s\t%s\t%s\t%0.3f\t%0.3f\t-\t-" % (name, alg1, alg2, a1_fmax.median(), a2_fmax.median())
            pval = 1
        sp_pval[s] = pval
        out_str += line+'\n'

    if kwargs['forceplot'] or not os.path.isfile(out_file):
        print("writing to %s" % (out_file))
        with open(out_file, 'w') as f:
            f.write(out_str)
    else:
        print("Would've written to %s. Use --forceplot to overwrite" % (out_file))
    print(out_str)

    # now check how many are significant
    num_sig = 0
    for s, pval in sp_pval.items():
        if pval*len(curr_species) < 0.05:
            num_sig += 1
    print("\t%d species with pval*%d < 0.05" % (num_sig, len(curr_species)))
    return sp_pval


def plot_fmax_eval(
    df_curr, out_file, measure='fmax', 
    sort_taxon_by_fmax=None, ann_stats=None, **kwargs):
    """
    *ann_stats*: Set of annotation types (either 'EXP', 'COMP' or 'IEA') for which a boxplot will be added to the right # of annotations
    *for_pub*: If true, the title at the top will not be included
    """
    if kwargs['forceplot'] or not os.path.isfile(out_file):
        print("Writing figure to %s" % (out_file))
    else:
        print("File would be written to %s" % (out_file))
        return
    algorithms = kwargs['algorithm']
    # add the species name to the boxplot
    species_labels = {}
    for s in df_curr['#taxon']:
        # the latex isn't working
        species_labels[s] = "%s. (%d)" % (
            f_settings.NAME_TO_SHORTNAME.get(selected_species[str(s)], str(s)),  #f_settings.NAME_TO_SHORTNAME2[selected_species[str(s)]],
            df_curr[df_curr['#taxon'] == s]['goid'].nunique())
    species_labels = pd.Series(species_labels)
    df_curr['species'] = df_curr['#taxon'].apply(lambda x: species_labels[x])
    df_curr = df_curr.sort_values(by=['species', 'Algorithm'], ascending=[True, False])
    # sort the species by fmax median of the first algorithm (sinksource)
    if sort_taxon_by_fmax is None:
        sort_taxon_by_fmax = df_curr[df_curr['Algorithm'] == kwargs['alg_names'][algorithms[0]]].groupby('#taxon')['fmax'].median().sort_values(ascending=False).index
        # now get the species order from the taxon order
    df_curr['#taxon'] = df_curr['#taxon'].astype("category")
    df_curr['#taxon'].cat.set_categories(sort_taxon_by_fmax, inplace=True)
    sort_by_med_fmax = df_curr.sort_values(by='#taxon')['species'].unique()

    xlabel = measure
    ylabel = 'Species (# GO Terms)'
#         fig, ax = plt.subplots(figsize=(6,10))
#         fig, ax = plt.subplots(figsize=(5,6))
    #if ann_stats is not None and ('COMP' in ann_stats or 'ELEC' in ann_stats):
    #    # make the figure taller to fit all the species
    #    fig, ax = plt.subplots(figsize=(3.5,8))
    #else:
    fig, ax = plt.subplots(figsize=(3.5,6))
    sns.boxplot(x=measure, y='species', order=sort_by_med_fmax, 
#                     hue='Algorithm', data=df_curr, orient='h')
               hue='Algorithm', hue_order=[kwargs['alg_names'][a] for a in algorithms], 
                data=df_curr, orient='h', fliersize=1.5,
               palette=plt_utils.my_palette)
    plt.title(kwargs['title'])
    plt.xlabel(xlabel, fontsize=12, weight="bold")
    plt.ylabel(ylabel, fontsize=12, weight="bold")
    # didn't work
#         locs, labels = plt.yticks()
#         plt.yticks(locs, [r'%s'%l for l in sort_by_med_fmax])
#         plt.legend(bbox_to_anchor=(1.2, 1.2))
    plt.legend(bbox_to_anchor=(.45, 1.0))
    ticks = np.arange(0,1.01,.1)
    plt.setp(ax, xticks=ticks, xticklabels=["%0.1f"%x for x in ticks])
    if ann_stats is not None:
        # add another plot to the right that is a bar plot of the # of non-iea (or exp) annotations
        right_ax = plt.axes([.95, .13, .2, .75])
        df_counts = pd.DataFrame()
        for ev_code in ann_stats:
            df_counts = pd.DataFrame([ann_stats[ev_code]]).T
            df_counts.columns = [ev_code]
        df_counts['#taxon'] = df_counts.index
        df_counts.index = range(len(df_counts.index))
        df_counts = pd.melt(df_counts, id_vars="#taxon", var_name="Ev. codes", value_name="# ann")
        df_counts['species'] = df_counts['#taxon'].apply(lambda x: species_labels[x])
        sns.barplot(y='species', x='# ann', hue="Ev. codes", order=list(sort_by_med_fmax), 
                    data=df_counts, orient='h', ax=right_ax, #palette="OrRd_d") 
                    color="#9c6d57") # use a brown color
        change_width(right_ax, 0.45)
        right_ax.set_ylabel("")
        right_ax.set_xlabel("# Ann")
        right_ax.set_yticks([])

    plt.savefig(out_file, bbox_inches='tight')
    plt.show()
    plt.close()
    return

# also plot a histogram of all the Fmax values
# cutoff1, cutoff2 = cutoffs[0]
# df_curr = cutoffs_data[0][0]
def plot_overview(cutoffs_data, cutoffs, out_dir, measure='fmax'):
    for i, (cutoff1, cutoff2) in enumerate(cutoffs):
        df_curr, _, _, _ = cutoffs_data[i]
        use_neg_str = "using negative examples" if with_neg != "" else ""
        out_pref = "%s/%d-%d-%s-%s%s-%s-%sa%s%s" % (out_dir, cutoff1, cutoff2, h, len(algorithms), with_neg, measure, unweighted, ss_alpha, maxi)
        title = "Recovery of %s ev. codes from %s ev. codes  %s\n for 19 pathogenic bacteria. - %s \n %d GO terms with %d <= # annotations < %d" % (
            eval_ev_codes, ev_codes, use_neg_str, version, df_curr['goid'].nunique(), cutoff1, cutoff2)
        print(df_curr['Algorithm'].unique())
#         print(df_curr.head())
#         print(title)
        overview_plots(df_curr, measure, title=title, out_pref=out_pref)

def overview_plots(df, measure, title="", xlabel=r'F$_{\mathrm{max}}$', 
                   out_pref=None, ax=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(5,4))
    # also plot a box and whisker plot of the algorithms
#     sns.factorplot(x=measure, y='Algorithm', data=df, ax=ax,
#                fliersize=1.5, order=[alg_name[a] for a in algorithms],
#                kind="box")
    sns.boxplot(x=measure, y='Algorithm', data=df, ax=ax,
               fliersize=1.5, order=[kwargs['alg_names'][a] for a in algorithms],
               palette=plt_utils.my_palette)
    ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
    ax.set_ylabel("")
    ax.set_title(title, fontsize=14, fontweight='bold')    
    ax.set_xlim(-0.025, 1.025)
    ticks = np.arange(0,1.01,.1)
    # ax.set_xticks(ticks, ["%0.1f"%x for x in ticks])
    plt.setp(ax, xticks=ticks, xticklabels=["%0.1f"%x for x in ticks])
    # add the intermediate xticks
#     locs, labels = plt.xticks()
    if out_pref is not None:
        out_file = out_pref + '-boxplot.pdf'
        print("Writing %s" % (out_file))
        plt.savefig(out_file, bbox_inches='tight')
#     plt.show()
#     plt.close()


# change the width of the extra bar plot to the side
def change_width(ax, new_value) :
    for patch in ax.patches :
        current_width = patch.get_height()
#         print(current_width)
        diff = current_width - new_value

        # we change the bar width
        patch.set_height(new_value)

        # we recenter the bar
        patch.set_y(patch.get_y() + diff * .5 + 0.08)
    
    
def get_results(versions, **kwargs):
    global cutoffs
    ev_code_results = {}
    for version in versions:
        for ev_codes, eval_ev_codes, string, h in itertools.product(
                *[kwargs['ev_codes'], kwargs['eval_ev_codes'], kwargs['string_nets']+[''], kwargs['hierarchy']]):
            exp_name = "%s-%s%s%s-%s%s" % (ev_codes, kwargs['size_range'], eval_ev_codes, '-%s'%string if string != '' else '', h.lower(), kwargs['pos_neg_str'])

            if 'string' not in version and \
               ('core' in exp_name or 'all' in exp_name or 'nontransferred' in exp_name):
                continue
            elif 'string' in version and \
                 ('core' not in exp_name and 'all' not in exp_name and 'nontransferred' not in exp_name):
                continue

            df_all = plt_utils.get_results_alg_params(version, exp_name, **kwargs)
            print(df_all.head(2))

            # TODO make an option for this
            cutoffs = [(10,1000)]
            if 's200' in version:
                cutoffs = [(10,5000)]
            cutoffs_data = split_results_by_cutoffs(df_all, cutoffs, exp_goid_prots, comp_goid_prots, iea_goid_prots)

            # 'recov' stands for recover 
            eval_str = "-recov-"+eval_ev_codes if len(eval_ev_codes) > 1 else eval_ev_codes
            out_dir = "outputs/viz/eval-loso/%s%s/%s" % (ev_codes, eval_str, version)
            utils.checkDir(out_dir)

            ev_code_results[(version, ev_codes, eval_ev_codes, h)] = (cutoffs_data, out_dir)
    return ev_code_results


def results_overview(ev_code_results, measures=['fmax'], **kwargs):
    if 'main_alg_opts' in kwargs:
        print("main_alg_opts:")
        for alg in kwargs['algorithm']:
            print("\t%s: %s" % (alg, kwargs['main_alg_opts'][alg]))
    for version, ev_codes, eval_ev_codes, h in sorted(ev_code_results, key=lambda x: (x[1], x[2], x[0])):
        cutoffs_data, out_dir = ev_code_results[(version, ev_codes, eval_ev_codes, h)]
        for df_curr, _,_,_ in cutoffs_data:
            num_vals = set()
            print(version, ev_codes, eval_ev_codes, h)
            # limit the goterms to those that are also present for SinkSource(?)
            for measure in measures:
                for alg in sorted(df_curr['Algorithm'].unique()):
                    vals = df_curr[df_curr['Algorithm'] == alg][measure]
                    print("\t%0.3f\t%s\t(%d sp-goterm pairs)" % (vals.median(), alg, len(vals)))
                    num_vals.add(len(vals))
            if len(num_vals) > 1:
                print("WARNING: not all algorithms have the same # of sp-goterm pairs")
                time.sleep(3)


def split_results_by_cutoffs(df_all, cutoffs, 
                             exp_goid_prots, comp_goid_prots, iea_goid_prots, 
                             overall_goid_prots=None):
    # now limit it to the current GO term split
    cutoffs_data = []
    for cutoff1, cutoff2 in cutoffs:
        print(cutoff1, cutoff2)
#         if split_by_overall_counts is True:
            # group the goids by the # of ann in all 19 species
            # this would be the # of annotations in the evaluation group
            # we want the # of annotations with the training evidence codes
            #num_ann = df_all[df_all['Algorithm'] == 'SinkSource'][['goid', '# ann']].groupby("goid").sum()
        if overall_goid_prots is not None:
            print("\tusing splits over all species")
            curr_goids = set([g for g, prots in overall_goid_prots.items() if len(prots) >= cutoff1 and len(prots) < cutoff2])
            df_curr = df_all[df_all['goid'].isin(curr_goids)]
            # also only keep those that have at least 10 in the evaluated species
            df_curr = df_curr[(df_curr['# ann'] >= 10)]
        else:
            print("\tusing splits for each individual species")
            col = '# test ann' if '# test ann' in df_all.columns else '# ann'
            df_curr = df_all[(df_all[col] >= cutoff1) & (df_all[col] < cutoff2)] 
        species_exp_counts = get_species_prot_counts(exp_goid_prots, species_to_uniprot, 
                                                     selected_species=df_curr['#taxon'].unique(), 
                                                     goids=set(df_curr['goid'].values))
        species_comp_counts = get_species_prot_counts(comp_goid_prots, species_to_uniprot, 
                                                      selected_species=df_curr['#taxon'].unique(), 
                                                      goids=set(df_curr['goid'].values))
        species_iea_counts = get_species_prot_counts(iea_goid_prots, species_to_uniprot, 
                                                     selected_species=df_curr['#taxon'].unique(), 
                                                     goids=set(df_curr['goid'].values))
        cutoffs_data.append((df_curr, species_exp_counts, species_comp_counts, species_iea_counts))
    return cutoffs_data


# also read in the prots annotated to each GO term
def parse_pos_neg_file(pos_neg_file, goterms=None):
    print("Reading positive examples for each goterm from %s" % (pos_neg_file))
    goid_prots = {}
    with open(pos_neg_file, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            goid, pos_neg_assignment, prots = line.rstrip().split('\t')[:3]
            if goterms and goid not in goterms:
                continue
            if int(pos_neg_assignment) == 1:
                prots = set(prots.split(','))
                goid_prots[goid] = prots
    return goid_prots


def get_species_prot_counts(goid_prots, species_to_uniprot, selected_species=None, goids=None):
    print("Getting annotation counts per species")
    if selected_species is None:
        selected_species = species_to_uniprot.keys()
    if goids is None:
        goids = goid_prots.keys()
    species_prot_counts = {}
    for s in selected_species:
        species_prots = species_to_uniprot[str(int(s))]
        #print("species: %s, %d prots" % (s, len(species_prots)))
#         total_prots = set()
        total_ann = 0
        for goid in goids:
            if goid not in goid_prots:
                continue
            s_prots = goid_prots[goid] & species_prots
#             total_prots.update(s_prots)
            total_ann += len(s_prots)
        species_prot_counts[s] = total_ann
        
    return species_prot_counts


def get_goid_prots(ev_codes="expc", ann_cutoff=50):
    pos_neg_files = [
            "inputs/pos-neg/%s/pos-neg-bp-%d-list.tsv" % (ev_codes, ann_cutoff), 
            "inputs/pos-neg/%s/pos-neg-mf-%d-list.tsv" % (ev_codes, ann_cutoff),
    ]
    exp_goid_prots = {}
    for pos_neg_file in pos_neg_files:
        exp_goid_prots.update(parse_pos_neg_file(pos_neg_file))
    return exp_goid_prots


def scatterplot_fmax(df_curr, curr_ev_codes, out_file, alg1="sinksource", alg2="localplus", **kwargs):
    # figure out which GO terms have the biggest difference for all species
    # plot a scatter plot of the differences across all GO terms

    # for some reason I have them flipped
    alg1 = kwargs['alg_names'].get(alg1, alg1)
    alg2 = kwargs['alg_names'].get(alg2, alg2)
    print("\nComparing fmax values of %s and %s" % (alg2, alg1))
    df_blast = df_curr[df_curr['Algorithm'] == alg2]
    # key: taxon, goid tuple. Value: fmax
    blast_scores = dict(zip(zip(df_blast['#taxon'], df_blast['goid']), df_blast['fmax']))
    df_ss = df_curr[df_curr['Algorithm'] == alg1]
    ss_scores = dict(zip(zip(df_ss['#taxon'], df_ss['goid']), df_ss['fmax']))

    goid_diffs = {}
    for taxon, goid in df_curr[['#taxon', 'goid']].values:
        if (taxon, goid) not in ss_scores and goid not in blast_scores:
            if kwargs['verbose']:
                print("WARNING: %s not in both SS and Blast Avg." % goid)
            continue
        if (taxon, goid) not in ss_scores:
            if kwargs['verbose']:
                print("WARNING: %s not in SS" % goid)
            continue
        if (taxon, goid) not in blast_scores:
            if kwargs['verbose']:
                print("WARNING: %s not in Blast Avg." % goid)
            continue
        goid_diff = ss_scores[(taxon, goid)] - blast_scores[(taxon, goid)]
    #     goid_diff = (ss_scores[(taxon, goid)] - blast_scores[(taxon, goid)]) / float(blast_scores[(taxon, goid)])
        goid_diffs[(taxon, goid)] = goid_diff
    print("\t%d %s, %d %s, %d diffs" % (len(blast_scores), alg2, len(ss_scores), alg1, len(goid_diffs)))

    # also load the summary of the GO terms
    summary_file = "inputs/pos-neg/%s/pos-neg-10-summary-stats.tsv" % (curr_ev_codes)
    df_summary = pd.read_csv(summary_file, sep='\t')
    goid_names = dict(zip(df_summary['GO term'], df_summary['GO term name']))
    goid_num_anno = dict(zip(df_summary['GO term'], df_summary['# positive examples']))
    goid_taxon_num_ann = dict(zip(zip(df_blast['#taxon'], df_blast['goid']), df_blast['# test ann']))

    diff_col = '%s - %s F-max'%(alg1, alg2)
    df = pd.DataFrame({'%s F-max'%alg2: blast_scores, '%s F-max'%alg1: ss_scores, diff_col: goid_diffs})
    # print the # and % where SS is >, <, and = local+
    print("%s > %s: %d (%0.3f)"% (alg1, alg2, len(df[df[diff_col] > 0]), len(df[df[diff_col] > 0]) / float(len(df))))
    print("%s < %s: %d (%0.3f)"% (alg1, alg2, len(df[df[diff_col] < 0]), len(df[df[diff_col] < 0]) / float(len(df))))
    print("%s = %s: %d (%0.3f)"% (alg1, alg2, len(df[df[diff_col] == 0]), len(df[df[diff_col] == 0]) / float(len(df))))

    print("\nTop 5 difference in f-max:")
    print(''.join(["%s, %s\t%s\t%s\t%s\n" % (
        t, g, goid_names[g], goid_num_anno[g], goid_diffs[t, g]) for t, g in sorted(
        goid_diffs, key=goid_diffs.get, reverse=True)[:5]]))

    if kwargs['forceplot'] or not os.path.isfile(out_file):
        print("Writing to %s" % (out_file))
        with open(out_file, 'w') as out:
            out.write("#taxon\tgoid\tname\t# ann\ttaxon # ann\t%s\t%s\tdiff\n" % (alg2, alg1))
            out.write(''.join(["%s\n" % (
                '\t'.join(str(x) for x in [t, g, goid_names[g], goid_num_anno[g], goid_taxon_num_ann[(t,g)],
                blast_scores[(t,g)], ss_scores[(t,g)], goid_diffs[(t,g)]])
            ) for t, g in sorted(goid_diffs, key=goid_diffs.get, reverse=True)]))
    else:
        print("Already exists: %s Use --forceplot to overwrite" % (out_file))

    # I can't get only the top histogram to change and not the right, so I have to make two copies of the file.
    # One with bins of 10 for the top, and the other with the regular # of bins
    for bins in [None, 10]:
        grid = sns.jointplot(x='%s F-max' % (alg1), y='%s - %s F-max'%(alg1, alg2), data=df,
                    stat_func=None, joint_kws={"s": 20}, marginal_kws=dict(bins=bins) if bins is not None else None,
                    )
        # plt.suptitle('%s, %s, %s-%s \n %d species' % (version, ev_codes, cutoff1, cutoff2, df_h_infl['#taxon'].nunique()))
        # plt.tight_layout()
        grid.fig.set_figwidth(5)
        grid.fig.set_figheight(5)
        out_file2 = out_file.replace('.tsv', '.pdf')
        if bins is not None:
            out_file2 = out_file.replace('.tsv', '-10bins.pdf')

        if kwargs['forceplot'] or not os.path.isfile(out_file2):
            print("Writing to %s" % (out_file2))
            plt.savefig(out_file2)
        else:
            print("Already exists: %s Use --forceplot to overwrite" % (out_file2))
        #plt.show()
        plt.close()


def parse_args(args):
    parser = plt_utils.setup_optparse()
    plt_utils.add_plot_opts(parser)

    # also add a couple of options specific to this script
    group = OptionGroup(parser, 'Stat-sig options')
    group.add_option('', '--alg1', default="sinksource",
                     help="Will test if/how much alg1 > alg2. --alg1 and --alg2 can also alternate opts (for example: 'sinksource a1.0-maxi1000'). Default: sinksource")
    group.add_option('', '--alg2', default="localplus",
                     help="Baseline to compare against. Default: localplus")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args(args)
    kwargs = plt_utils.validate_opts(opts) 

    return kwargs
    

if __name__ == "__main__":
    kwargs = parse_args(sys.argv)
    versions = kwargs.pop('version')

    # variables used to get the # of annotations for each ev code type.
    # TODO organize these better 
    exp_goid_prots = get_goid_prots(ev_codes="expc", ann_cutoff=10) 
    comp_goid_prots = get_goid_prots(ev_codes="comp", ann_cutoff=5) 
    iea_goid_prots = get_goid_prots(ev_codes="iea", ann_cutoff=10) 

    main(versions, **kwargs)
