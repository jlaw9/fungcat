#!/usr/bin/python

print("Importing libraries")

import os
import sys
import glob
import shutil
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


# use the paired pallet to pair the ovn and ova
current_palette = sns.color_palette("Paired")
sns.set_palette(current_palette)


COLORS = {
        "sinksource-ovn": current_palette[0],
        "sinksource-ova": current_palette[1],
        "local-ovn": current_palette[2],
        "local-ova": current_palette[3],
        "fun-flow-ovn": current_palette[4],
        "genemania-ova": current_palette[5],
}

SHAPES = {
        "local-ova": 'o',
        "local-ovn": 'X',
        "fun-flow-ovn": 's',
        "sinksource-ova": '>',
        "sinksource-ovn": 'v',
        "genemania-ova": '^',
}

VERSION_SHAPES = ['o', 'X','s','>','v','^']

NAMES = {
        "local-ova": "Blast",
        "local-ovn": 'Blast+',
        "fun-flow-ovn": 'Functional Flow',
        "sinksource-ova": 'SinkSource',
        "sinksource-ovn": 'SinkSource+',
        "genemania-ova": 'GeneMANIA',
}

# Guide to GO evidence codes: http://geneontology.org/page/guide-go-evidence-codes
ALL_EVIDENCE_CODES = ['EXP','IDA','IPI','IMP','IGI','IEP','ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA','TAS','NAS','IC','ND','IEA']
EXP_EVIDENCE_CODES = ['EXP','IDA','IPI','IMP','IGI','IEP']
COMP_EVIDENCE_CODES = ['ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA']
AUTH_EVIDENCE_CODES = ['TAS','NAS']
CUR_EVIDENCE_CODES = ['IC']


def main(versions, exp_names, algorithms, selected_terms, depth_by_numanno=False, pdf=False, forced=False):
    """
    *versions*: versions will be compared to each other on the same plot (fmeasure)
    *selected_terms*: will create one prec-rec plot for each goid. Must be only the numbers (i.e., no leading GO:00)
    """

    version_cv_results, goids_info = get_cv_results(versions, exp_names, algorithms, selected_terms)

    version_prec_rec_results = {version: defaultdict(dict) for version in versions}
    # quick hack to get the pr rec results out of the CV results 
    # used by many other functions
    for version in version_cv_results:
        for goid in version_cv_results[version]:
            for alg in version_cv_results[version][goid]:
                pr_rec = [(pr, rec) for pr, rec, fp in version_cv_results[version][goid][alg]]
                version_prec_rec_results[version][goid][alg] = pr_rec

    # organize the figures by the exp names being combined
    out_dir = "outputs/viz/eval/%s" % ('-'.join(exp_names))
    utils.checkDir(out_dir)

    # make a scatter plot 
    if len(versions) == 2:
    #if False:
        out_file = "%s/fmax-%s-%s-numanno.png" % (out_dir, versions[0], versions[1])
        # TODO this is handled elsewhere now
        if (forced is False and os.path.isfile(out_file)) and False:
            print( "%s already exists. Used --forced to overwrite it" % (out_file))
        else:
            title = "Scatterplot of SinkSource F-max results\ncolored by # of genes annotated"
            scatterplot_versions_with_category(version_prec_rec_results, goids_info, algorithm="sinksource-ova",
                                               title=title, out_file=out_file, pdf=pdf)
            # TODO make this an option or something
            out_file = "%s/fmax-%s-%s.png" % (out_dir, versions[0], "2018_02-stacking-rf")
            goterm_fmax_file = "/data/jeff-law/projects/fungcat-function-prediction/inputs/datasink_rf_fmax.csv"
            title = "Scatterplot of F-max values for SinkSource vs Datasink RF \ncolored by # genes annotated"
            scatterplot_other(version_prec_rec_results, goterm_fmax_file, goids_info,
                              version="2017_10-seq-sim", title=title, out_file=out_file, pdf=pdf)

    title = "Scatterplot of SinkSource F-max results\ncolored by ratio of exp. to comp. evidencecodes"
    out_file = "%s/fmax-%s-%s-delta-evidencecodes.png" % (out_dir, versions[0], versions[1])
    scatterplot_versions_with_evidencecodes(
            version_prec_rec_results, goids_info, algorithm="sinksource-ova",
            title=title, out_file=out_file, pdf=pdf, forced=forced)

    # TODO add an option to differentiate these different analyses
    if len(exp_names) > 1:
        # quit here because the plots below add the experiment name to the plot 
        print("Quitting")
        return
    exp_name = exp_names[0]

    if not opts.goterm_plots:
        print("Skipping creating fmax and prec/rec plots for individual GO terms")
        return

    # plot the f-max with versions on the same plot first
    for selected_term in goids_info:
        goid_info = goids_info[selected_term]
        #print(selected_term)
        #print(goid_info)
        goname = goid_info['name'] if len(goid_info['name']) < 40 else "%s..."%(goid_info['name'][:40])
        title = "Function '%s' (%s)\ncategory: %s, depth: %s, genes: %s\nexperiment: %s" \
                % (goname, selected_term, goid_info['category'], goid_info['depth'], goid_info['genes'], exp_name) 
        out_file_name = "%s-fmax.png" % (goid_info['name'].replace(' ', '-'))
        out_file = "%s/%s" % (out_dir, out_file_name)
        goid_pr = {}
        for version in versions:
            goid_pr[version] = version_prec_rec_results[version][selected_term]

        plot_fmax(goid_pr, algorithms, title, out_file=out_file, pdf=pdf, forced=forced)

    # now plot the prec/rec
    for version in versions:
        INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains = f_settings.set_version(version)
        out_dir = "%s/all/viz/eval/%s" % (RESULTSPREFIX, exp_name)
        utils.checkDir(out_dir)
        for selected_term in goids_info:
            goid_info = goids_info[selected_term]
            # also plot the precision and recall
            out_file_name = "%s-pr.png" % (goid_info['name'].replace(' ', '-'))
            out_file = "%s/%s" % (out_dir, out_file_name)
            goname = goid_info['name'] if len(goid_info['name']) < 40 else "%s..."%(goid_info['name'][:40])
            title = "Function '%s' (%s)\ncategory: %s, depth: %s, genes: %s\nversion: %s, experiment: %s" \
                    % (goname, selected_term, goid_info['category'], goid_info['depth'], goid_info['genes'], version, exp_name) 
            # now plot the prec/rec
            plot_prec_rec(version_prec_rec_results[version][selected_term], title, 
                    out_file=out_file, pdf=pdf, forced=forced)
            ## now copy the plots to a more friendly viz directory
            #out_dir2 = "outputs/viz/%s/%s" % (exp_name, version)
            #print "Copying plot to %s" % (out_dir2)
            #utils.checkDir(out_dir2)
            #shutil.copy(out_file, out_dir2)

    if depth_by_numanno:
        # plot the f-max with versions on the same plot first
        out_dir = "outputs/stats"
        utils.checkDir(out_dir)
        out_file = "%s/%s-depth-by-numanno.png" % (out_dir, version)
        title = "Depth of GO term in GO DAG by # of Annotations\nversion: %s" % (version)

        plot_depth_by_numanno(goids_info, title, out_file, pdf=pdf)


def get_annotation_evidencecodes(gain_file):
    print("Reading annotations from GAIN file %s. Assuming annotations have already been propogated up the GO DAG" % (gain_file))

    # dictionary with key: goterm ID, val: set of proteins annotated to the goterm ID 
    evidencecodes_per_goid = {}

    # columns: 'orf', 'goid', 'hierarchy', 'evidencecode', 'annotation type' (described here: http://bioinformatics.cs.vt.edu/~murali/software/biorithm/gain.html)
    df = pd.read_csv(gain_file, sep='\t')
    print("%d proteins have at least 1 annotation" % (len(set(df['orf'].tolist()))))
    #print("%d goids" % (len(set(df['goid'].tolist()))))
    # change the hierarchy to be upper case
    df['hierarchy'] = df['hierarchy'].apply(lambda x: x.upper())
    # also add the leading GO:000 to the ID
    df['goid'] = df['goid'].apply(lambda x: "GO:" + "0"*(7-len(str(x))) + str(x))
    #goid_per_prot_by_h = {}
    #for h in ["C", "F", "P"]:
        #print(h)
        #df_h = df[df["hierarchy"] == h]
        #goid_per_prot_by_h[h] = {orf: set(goids['goid'].tolist()) for orf, goids in df_h[['orf', 'goid']].groupby('orf')}
        #print(h, len(goid_per_prot_by_h[h]))
    evidencecodes_per_goid = {goid: evidencecodes['evidencecode'].tolist() for goid, evidencecodes in df[['evidencecode', 'goid']].groupby('goid')}
    #print(len(prot_per_goid))
    #goid_to_hierarchy = dict(zip(df['goid'], df['hierarchy']))
    #print(len(goid_to_hierarchy))

    return evidencecodes_per_goid


def get_cv_results(versions, exp_names, algorithms, selected_terms):
    """

    returns a 3 level dictionary with a list of precision, recall, fp_rate tuples
      version:
        goid:
          algorithm:
            [ (prec, rec, fp_rate), ... ]
    """

    print("Parsing CV results for %d versions, %d exp_names, %d algorithms, and %d selected_terms" % (len(versions), len(exp_names), len(algorithms), len(selected_terms)))
    version_cv_results = {version: defaultdict(dict) for version in versions}
    # GOID info does not need to be per version because the information will be the same
    all_goids_info = {}
    for version in versions:
        for exp_name in exp_names:
            INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains = f_settings.set_version(version)

            # get the cross-validation file(s)
            gain_results_dir = "%s/all" % (RESULTSPREFIX)
            #print("Plotting CV results for %d terms: %s" % (len(selected_terms), str(selected_terms)))
            #for selected_term in selected_terms:
                #gain_results_combined = "%s/db-%s-cv.txt" % (out_dir, exp_name)
                #print "Combining all of the different algorithm's results to %s, and plotting the precision by recall" % (gain_results_combined)
                #with open(gain_results_combined, 'w') as out:
                # load the prediction results of each algorithm
            alg_goid_cv_results, goids_info = parse_cv_file(gain_results_dir, algorithms, exp_name, selected_terms)
                
            for alg in algorithms:
                for goid in selected_terms:
                    if goid not in goids_info:
                        continue
                    # prepend GO: with leading 0s
                    goid_fixed = "GO:" + "0"*(7-len(str(goid))) + str(goid)
                    version_cv_results[version][goid_fixed][alg] = alg_goid_cv_results[alg][goid]
                    #goids_info[goid] = goid_info
                    all_goids_info[goid_fixed] = goids_info[goid]

    return version_cv_results, all_goids_info 


def parse_cv_file(gain_results_dir, algorithms, exp_name, selected_terms):
    # each algorithm contains a dictionary of go term: cv results (list of prec, rec, fpr tuples)
    alg_goid_cv_results = {alg: defaultdict(list) for alg in algorithms}
    # its the same for each algorithm
    #goid_info = {alg: {} for alg in algorithms}
    goid_info = {}
    found_terms = set()
    for alg in algorithms:
        #gain_results_file = "%s/%s/pred/db-%s-cv.txt" % (gain_results_dir, alg, exp_name)
        gain_results_file = "%s/%s/%s/db-%s-cv.txt" % (gain_results_dir, alg, exp_name, exp_name)
        print( "\treading %s for %d terms" % (gain_results_file, len(selected_terms)))
        with open(gain_results_file, 'r') as f:
            # the cv file always has two blank lines
            line = f.readline()
            line = f.readline()
            line = f.readline()
            while line != "":
                #print line
                # header line 1
                info_line = line.rstrip().split('  ')[-1]
                info = {}
                for pair in info_line.split(', '):
                    try:
                        #name,val = pair.split(' ')
                        name = pair[:pair.index(' ')]
                        val = pair[pair.index(' ')+1:]
                        #print 'name:', name, 'val:', val
                        info[name] = val
                    except ValueError:
                        print("Failed to parse pair %s. Skipping. Info line:" % (str(pair)))
                        print(info_line)
                        continue
                    
                #info = {name:val for pair in info_line.split(', ') for name,val in pair.split(' ')}
                # genes isn't right...
                info['genes'] = info_line.split(', ')[-1].split(' ')[0]
                goid = info['id']
                if goid not in selected_terms:
                    #print "skipping %s" % (goid)
                    # read lines until the next 2 blank lines
                    while line != "":
                        line = f.readline().rstrip()
                else:
                    found_terms.add(goid)
                    #print("\t%s" % (line.rstrip()))
                    #print("\treading %s" % (goid))
                    goid_info[goid] = info
                    header = f.readline()
                    line = f.readline().rstrip()
                    while line != "":
                        #print line
                        line = line.split('\t')
                        #Confidence Cutoff  Desired Recall  Actual Recall   Precision   FP Rate TP  FP  TN  FN
                        rec = float(line[2])
                        prec = float(line[3])
                        fp_rate = float(line[4])
                        alg_goid_cv_results[alg][goid].append((prec, rec, fp_rate)) 
                        line = f.readline().rstrip()

                line = f.readline()
                line = f.readline()

    if len(found_terms) != len(selected_terms):
        print( "\tWARNING: Found data for only %d out of %d GO terms" % (len(found_terms), len(selected_terms)))

    return alg_goid_cv_results, goid_info


def plot_depth_by_numanno(goids_info, title, out_file="test.png", pdf=False):
    goids = goids_info.keys()
    depths = [int(goids_info[goid]['depth']) for goid in goids]
    num_anno = [int(goids_info[goid]['genes']) for goid in goids]

    fig, ax = plt.subplots(figsize=(6,5))

    plt.scatter(depths, num_anno, color=current_palette[1])

    plt.xlabel("Depth of term in GO DAG")
    plt.ylabel("# Annotations")
    plt.title(title)

    plt.xticks(range(min(depths),max(depths)+1))
    ax.set_xticklabels(range(min(depths),max(depths)+1))

    savefig(fig, out_file, pdf=pdf) 


def scatterplot_other(
        version_prec_rec_results, goid_fmax_file, goids_info,
        version="", algorithm="sinksource-ova", title="Scatterplot of F-max results for SinkSource",
        out_file='src/test.png', pdf=False):

    # get the fmax of the given algorithm for each version
#    my_goid_fmax = defaultdict(dict)
#    for goid in prec_rec_results:
#        alg_prec_rec = prec_rec_results[goid][algorithm]
#        f_measures = []
#        for p,r in alg_prec_rec:
#            harmonic_mean = (2*p*r)/(p+r)
#            f_measures.append(harmonic_mean)
#        goid_fmax[version][goid] = max(f_measures)
    version_goid_fmax = defaultdict(dict)
    for goid in version_prec_rec_results[version]:
        alg_prec_rec = version_prec_rec_results[version][goid][algorithm]
        f_measures = []
        for p,r in alg_prec_rec:
            harmonic_mean = (2*p*r)/(p+r)
            f_measures.append(harmonic_mean)
        version_goid_fmax[version][goid] = max(f_measures)
    other_goid_fmax = utils.readDict(goid_fmax_file, sep=',')
    #print(other_goid_fmax)
    other_goid_fmax = {goid: float(fmax) for goid, fmax in other_goid_fmax.items()}
    #version_prec_rec_results['2017_01-stacking']
    other_version = "2018_02-stacking-random-forest"

    df = pd.DataFrame({version: version_goid_fmax[version], other_version: other_goid_fmax})
    category_names = {"f": "MF", "p": "BP"}
    goid_category = {goid: category_names[goids_info[goid]['category']] for goid in goids_info} 
    #print(goid_category)
    df['category'] = pd.Series(goid_category) 
    df['depth'] = pd.Series({goid: goids_info[goid]['depth'] for goid in goids_info})
    df['name'] = pd.Series({goid: goids_info[goid]['name'] for goid in goids_info})
    df['num_annotations'] = pd.Series({goid: goids_info[goid]['genes'] for goid in goids_info})
    # remove rows that don't have values for both
    df.dropna(how='any',inplace=True)
    df.to_csv(out_file.replace('.png', '.csv'))
    #print(df)
    #sys.exit()

    # now make the scatterplot
    fig, ax = plt.subplots(figsize=(6,5))
    #v1, v2 = version_goid_fmax.keys()
    #sns.pairplot(x_vars=[version], y_vars=[other_version], data=df, hue='Category', size=5, palette=sns.color_palette("Set1"))
    cmap = sns.cubehelix_palette(as_cmap=True)
    points = ax.scatter(df[version], df[other_version], c=df['num_annotations'], cmap=cmap)
    fig.colorbar(points)
    plt.plot([0,1])

    plt.ylim(ymax=1.00, ymin=0)
    plt.xlim(xmax=1.00, xmin=0)

    plt.xlabel(version)
    plt.ylabel(other_version)
    plt.title(title)

    savefig(fig, out_file, pdf=pdf) 

    
def scatterplot_versions_with_category(
        version_prec_rec_results, goids_info,
        algorithm="sinksource-ova", title="Scatterplot of F-max results for SinkSource",
        out_file='src/test.png', pdf=False):
    """
    """
    # get the fmax of the given algorithm for each version
    version_goid_fmax = defaultdict(dict)
    for version in version_prec_rec_results:
        for goid in version_prec_rec_results[version]:
            alg_prec_rec = version_prec_rec_results[version][goid][algorithm]
            #print(version, goid, algorithm, alg_prec_rec)
            f_measures = []
            for p,r in alg_prec_rec:
                harmonic_mean = (2*p*r)/(p+r)
                f_measures.append(harmonic_mean)
            version_goid_fmax[version][goid] = max(f_measures)

    df = pd.DataFrame(version_goid_fmax)
    category_names = {"f": "MF", "p": "BP"}
    goid_category = {goid: category_names[goids_info[goid]['category']] for goid in goids_info} 
    #print(goid_category)
    df['Category'] = pd.Series(goid_category) 
    df['depth'] = pd.Series({goid: goids_info[goid]['depth'] for goid in goids_info})
    df['name'] = pd.Series({goid: goids_info[goid]['name'] for goid in goids_info})
    df['num_annotations'] = pd.Series({goid: goids_info[goid]['genes'] for goid in goids_info})
    df.to_csv(out_file.replace('.png', '.csv'))
    #print(df)
    #sys.exit()

    # now make the scatterplot
    fig, ax = plt.subplots(figsize=(6,5))
    v1, v2 = version_goid_fmax.keys()
    #sns.pairplot(x_vars=[v2], y_vars=[v1], data=df, hue='Category', size=5, palette=sns.color_palette("Set1"))
    cmap = sns.cubehelix_palette(as_cmap=True)
    points = ax.scatter(df[v2], df[v1], c=df['num_annotations'], cmap=cmap)
    fig.colorbar(points)
    plt.plot([0,1])
    #colors = {}
    #colors[df['category']=="BP"] = current_palette[0]
    #colors[df['category']=="MF"] = current_palette[2]
    #print colors
    #df.plot.scatter(x=df[v1], y=df[v2], c=colors)
    #category_colors = {"f": current_palette[0], "p": current_palette[1]}
    #for goid in goids_info:
    #    color = category_colors[goids_info[goid]['category']]
    #    label = category_names[goids_info[goid]['category']]
    #    ax.plot(version_goid_fmax[v1][goid], version_goid_fmax[v2][goid], ms=12, color=color, label=label)
#
    plt.ylim(ymax=1.00, ymin=0)
    plt.xlim(xmax=1.00, xmin=0)

    plt.xlabel(v2)
    plt.ylabel(v1)
    plt.title(title)

    savefig(fig, out_file, pdf=pdf) 


def scatterplot_versions_with_evidencecodes(
        version_prec_rec_results, goids_info,
        algorithm="sinksource-ova", title="Scatterplot of F-max results for SinkSource",
        out_file='src/test.png', pdf=False, forced=False):
    """
    color by the ratio of experimental to computational evidence codes
    """
    # get the fmax of the given algorithm for each version
    version_goid_fmax = defaultdict(dict)
    for version in version_prec_rec_results:
        for goid in version_prec_rec_results[version]:
            alg_prec_rec = version_prec_rec_results[version][goid][algorithm]
            #print(version, goid, algorithm, alg_prec_rec)
            f_measures = []
            for p,r in alg_prec_rec:
                harmonic_mean = (2*p*r)/(p+r)
                f_measures.append(harmonic_mean)
            version_goid_fmax[version][goid] = max(f_measures)

    gain_ann_file = f_settings.GOA_ALL_FUN_FILE_NOIEA
    evidencecodes_per_goid = get_annotation_evidencecodes(gain_ann_file)
    goid_ratio_exp_comp = {}
    for goid in goids_info:
        num_exp = len([ec for ec in evidencecodes_per_goid[goid] \
                if ec in EXP_EVIDENCE_CODES or ec in AUTH_EVIDENCE_CODES])
                #if ec in EXP_EVIDENCE_CODES])
        num_comp = len([ec for ec in evidencecodes_per_goid[goid] \
                if ec in COMP_EVIDENCE_CODES])
        if num_comp != 0:
            goid_ratio_exp_comp[goid] = float(num_exp) / float(num_comp) 
        else:
            print "0 computational evidence codes for %s" % (goid)
            goid_ratio_exp_comp[goid] = -1

    df = pd.DataFrame(version_goid_fmax)
    category_names = {"f": "MF", "p": "BP"}
    goid_category = {goid: category_names[goids_info[goid]['category']] for goid in goids_info} 
    #print(goid_category)
    df['Category'] = pd.Series(goid_category) 
    df['depth'] = pd.Series({goid: goids_info[goid]['depth'] for goid in goids_info})
    df['name'] = pd.Series({goid: goids_info[goid]['name'] for goid in goids_info})
    df['num_annotations'] = pd.Series({goid: goids_info[goid]['genes'] for goid in goids_info})
    df['ratio_exp_comp'] = pd.Series(goid_ratio_exp_comp)
    df.to_csv(out_file.replace('.png', '.csv'))
    #print(df)
    #sys.exit()

    # now make the scatterplot
    fig, ax = plt.subplots(figsize=(6,5))
    v1, v2 = version_goid_fmax.keys()
    #sns.pairplot(x_vars=[v2], y_vars=[v1], data=df, hue='Category', size=5, palette=sns.color_palette("Set1"))
    #cmap = sns.cubehelix_palette(as_cmap=True)
    points = ax.scatter(df['ratio_exp_comp'], df[v2] - df[v1], color=current_palette[1])
    #fig.colorbar(points)
    plt.axvline(1)
    #colors = {}
    #colors[df['category']=="BP"] = current_palette[0]
    #colors[df['category']=="MF"] = current_palette[2]
    #print colors
    #df.plot.scatter(x=df[v1], y=df[v2], c=colors)
    #category_colors = {"f": current_palette[0], "p": current_palette[1]}
    #for goid in goids_info:
    #    color = category_colors[goids_info[goid]['category']]
    #    label = category_names[goids_info[goid]['category']]
    #    ax.plot(version_goid_fmax[v1][goid], version_goid_fmax[v2][goid], ms=12, color=color, label=label)
#
    #plt.ylim(ymax=1.00, ymin=0)
    #plt.xlim(xmax=1.00, xmin=0)

    plt.xlabel("Ratio of Exp / Comp")
    plt.ylabel("%s fmax - %s fmax" % (v2, v1))
    plt.title(title)

    savefig(fig, out_file, pdf=pdf, forced=forced) 


def plot_fmax(version_prec_rec_results, algorithms, title, 
        out_file='src/test.png', pdf=False, forced=False):
    fig, ax = plt.subplots(figsize=(6,5))
    #alg_order = [
    #    "local-ova",
    #    "sinksource-ova",
    #    "genemania-ova",
    #    "local-ovn",
    #    "sinksource-ovn",
    #    "fun-flow-ovn",
    #    ]

    #version_alg_fmax = defaultdict(dict)
    alg_fmax = {}
    for i, version in enumerate(sorted(version_prec_rec_results)):
        x = 0
        algs = []
        prec_rec_results = version_prec_rec_results[version]
        for alg in algorithms:
            if alg not in prec_rec_results or len(prec_rec_results[alg]) == 0:
                continue
            algs.append(alg)
            pr = prec_rec_results[alg]
            f_measures = []
            for p,r in pr:
                harmonic_mean = (2*p*r)/(p+r)
                f_measures.append(harmonic_mean)
            alg_fmax[alg] = max(f_measures)

            x += 1
            # now plot them
            # make everything the same shape for a given version
            markerstyle = VERSION_SHAPES[i]
            markersize = 10
            color = COLORS[alg]

            ax.plot(x, alg_fmax[alg], markerstyle, ms=markersize, color=color)
    #df = pd.DataFrame([[NAMES[alg], alg_fmax[alg]] for alg in alg_order if alg in alg_fmax], columns=["Algorithm", "F-max"])
    #sns.stripplot(x="Algorithm", y="F-max", data=df, size=10)

    plt.xticks(range(1,len(algs)+1))
    ax.set_xticklabels([NAMES[alg] for alg in algs], rotation=10)
    plt.ylim(ymax=1.00, ymin=0)

    ax.set_ylabel('F-max', size=12)
    ax.set_title(title)

    # add a legend
    handles = []
    for i, version in enumerate(sorted(version_prec_rec_results)):
        shape = mlines.Line2D([], [], linestyle='None', marker=VERSION_SHAPES[i],
                                      markersize=12, color='gray', label=version)
        #patch = mpatches.Patch(shape=VERSION_SHAPES[i], label=version)
        handles.append(shape)
    legend = ax.legend(handles=handles, frameon=True)
    # Put a nicer background color on the legend.
    legend.get_frame().set_facecolor('#ffffff')
    legend.get_frame().set_linewidth(2.0)

    savefig(fig, out_file, pdf=pdf, forced=forced) 

    return

    
def plot_prec_rec(prec_rec_results, title, out_file='src/test.png', pdf=False, forced=False):
    #plt.gca().set_color_cycle(
    # use the paired pallet to pair the ovn and ova
    fig, ax = plt.subplots(figsize=(6,5))
    alg_order = [
        "local-ova",
        "sinksource-ova",
        "genemania-ova",
        "local-ovn",
        "sinksource-ovn",
        "fun-flow-ovn",
        ]

    for alg in alg_order:
        if alg not in prec_rec_results:
            continue
        pr = prec_rec_results[alg]
        # make everything a line as well as a shape for each prec-rec point
        markerstyle = SHAPES[alg] + '-'
        markersize = 10
        linewidth = 1
        color = COLORS[alg]
        label = NAMES[alg]
        #print alg, pr, markerstyle, color, label

        ax.plot([r for p,r in pr], [p for p,r in pr], 
                markerstyle, ms=markersize, lw=linewidth, 
                color=color, label=label)

    ax.set_xlabel('Recall', size=12)
    ax.set_ylabel('Precision', size=12)
    ax.set_title(title)

    legend = ax.legend(frameon=True)
    # Put a nicer background color on the legend.
    legend.get_frame().set_facecolor('#ffffff')
    legend.get_frame().set_linewidth(2.0)
    # if the title is too long, the plot does weird things
    plt.tight_layout()
    plt.ylim(ymax=1.00, ymin=0)

    savefig(fig, out_file, pdf=pdf, forced=forced) 

    return


def savefig(fig, out_file, pdf=False, forced=False):
    if forced is False and os.path.isfile(out_file):
        print("%s already exists. Used --forced to overwrite it" % (out_file))
        return
    print('\twriting to %s' % (out_file))
    plt.savefig(out_file)
    if pdf:
        out_file = out_file.replace('.png', '.pdf')
        print('\twriting to %s' % (out_file))
        plt.savefig(out_file, bbox_inches='tight')


# TODO modify the utils post_to_graphspace.py parseArgs so I can call it instead
def parseArgs(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--version',type='string', action='append',
                      help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    parser.add_option('', '--exp-name', type='string', action='append',
                      help="Experiment name to use when running GAIN. Can specify multiple experiment names to get results for")
    parser.add_option('-a', '--algorithm', action="append", 
                      help="Algorithm(s) to use when running GAIN. If not specified, many algorithms will be used. Options: '%s'" % ("', '".join(f_settings.ALGORITHM_OPTIONS.keys())))
    parser.add_option('-s', '--selected-term', type='string', action="append",
                      help="Selected goid to plot. Must be only the numbers (i.e., no leading GO:00). Can be specified multiple times")
    parser.add_option('', '--only-functions', type='string', action='append',
                      help="File containing GOIDs to plot")
    parser.add_option('', '--depth-by-numanno', action='store_true', default=False,
                      help="Plot the number of annotations per GO term by the depth of that term in the GO DAG")
    parser.add_option('', '--goterm-plots', action='store_true', default=False,
                      help="Create an fmax and prec/rec plot for each GO term")
    parser.add_option('', '--test-sig', action='store_true', default=False,
                      help="Run Kruskal-Wallis and Mann-Whitney U tests on boxplots (plot_goterm_sig.py only)")
    parser.add_option('', '--pdf', action='store_true', default=False,
                      help="Also plot a pdf")
    parser.add_option('', '--forced', action='store_true', default=False,
                      help="Force generating plots if they already exist")

    (opts, args) = parser.parse_args()

    if opts.exp_name is None:
        print("--exp-name required")
        sys.exit(1)

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)

    return opts, args

if __name__ == '__main__':
    opts, args = parseArgs(sys.argv)

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
    main(opts.version, opts.exp_name, algorithms, only_functions, depth_by_numanno=opts.depth_by_numanno, pdf=opts.pdf, forced=opts.forced)
