
from optparse import OptionParser, OptionGroup
from collections import defaultdict
import os
import sys 
from tqdm import tqdm
import itertools
sys.path.append('src')
import utils.file_utils as utils
import fungcat_settings as f_settings
import algorithms.gain_scipy.alg_utils as alg_utils
# import plotting libraries
import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import time
import json


# tried to be fancy :P
# colors: https://coolors.co/ef6e4a-0ec9aa-7c9299-5d88d3-96bd33
my_palette = ["#EF6E4A", "#0EC9AA", "#7C9299", "#5D88D3", "#96BD33", "#937860", "#efd2b8"]
# for comparing sinksource with local
#my_palette = ["#EF6E4A", sns.xkcd_rgb["deep sky blue"], "#96BD33", "#937860", "#efd2b8"]
#my_palette = [sns.xkcd_rgb["orange"], sns.xkcd_rgb["azure"], "#96BD33", "#937860", "#efd2b8"]

# these are the main options used to determine which parameters to use for the algorithm
# (unless the plot is specifically comparing parameters, then all will be used)
# for now, just check if the entire string is in the results file. If it is, then use it
# TODO use the JSON file
MAIN_WEIGHT_TYPE = {
    "birgrank": "swsn",
    "sinksource": "swsn",
    "genemania": "swsn",
    "localplus": "swsn",
    }
# TODO make an option for this
MAIN_ALG_OPTS = {
    'sinksource': '-l0-a0_95-eps0_0-maxi10.txt',  # TODO maxi10 at the end also matches maxi1000
    #'sinksource': '-l0-a0_95-eps0_0-maxi20.txt',  # TODO maxi10 at the end also matches maxi1000
    'genemania': "-l0-a0_95-eps0_0-maxi20-tol1e-05",  # TODO remove a, eps and maxi from the output file name
    'birgrank': "-l0-a0_9-eps0_0001-maxi1000-t0_5-m0_5-l0_01",
    'localplus': "-l0-a0_95-eps0_0-maxi10.txt",  # TODO remove a, eps and maxi from the output file name
    }
ALG_NAMES = {
    'localplus': 'Local', 'local': 'Local-',
    'sinksource': 'SinkSource', 'sinksourceplus': 'SinkSource+',
    'genemania': 'GeneMANIA',  
    'birgrank': 'BirgRank', 'aptrank': 'AptRank OneWay',
    'sinksource-squeeze': 'SinkSource Squeeze', 
    }


def main(versions, **kwargs):
    kwargs['main_weight_type'] = MAIN_WEIGHT_TYPE 
    kwargs['main_alg_opts'] = MAIN_ALG_OPTS 
    kwargs['alg_names'] = ALG_NAMES 

    ev_code_results = get_results(versions, **kwargs)
    results_overview(ev_code_results, **kwargs) 


def get_results(versions, **kwargs):
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

            df_all = get_results_alg_params(version, exp_name, **kwargs)
            #print(df_all.head())
            ev_code_results[(version, exp_name)] = df_all
    return ev_code_results


#def get_results(version, weight_type, **kwargs):
#main_alg_opts=None, main_weight_type=None, alg_names=None, 
def get_results_alg_params(version, exp_name, **kwargs):
    """ Get the fmax results for a given set of algorithm parameters
    """
    # see if these are set
    main_alg_opts = kwargs.pop('main_alg_opts', None)
    main_weight_type = kwargs.pop('main_weight_type', None)
    alg_names = kwargs.pop('alg_names', None)

    results_files = set()
    df_all = pd.DataFrame()
    for alg in kwargs['algorithm']:
        for weight_type, alpha, theta, mu, br_lambda, eps, maxi, tol in itertools.product(
                *[kwargs['weight_type'], kwargs['alpha'], kwargs['theta'], kwargs['mu'], kwargs['br_lambda'], kwargs['eps'], kwargs['max_iters'], kwargs['tol']]):

            # setup the weight_type
            if 'string' in version:
                if alg == 'birgrank' and weight_type in ['', 'gm2008']:
                    weight_type = 'swsn'
            # if the version doesn't have STRING in it, then don't use the weight str
            else:
                weight_type = "" 
            # add the dash to the weight str
            weight_str = "-%s"%weight_type if weight_type != "" else "" 

            results_file = alg_utils.get_filepath(
                version, alg, exp_name, weight_str=weight_str, exp_type='loso', 
                alpha=alpha, theta=theta, mu=mu, br_lambda=br_lambda, eps=eps, max_iters=maxi, tol=tol,
                sinksourceplus_lambda=kwargs['sinksourceplus_lambda']) 
            # make sure to not read files multiple times
            # for example: if alg is sinksource, but there are multiple 'tol' parameters, 
            # they will be ignored when getting the results file
            if results_file in results_files:
                continue
            results_files.add(results_file)

            if kwargs['verbose']:
                print("Reading results from %s" % (results_file))
            if not os.path.isfile(results_file):
                if kwargs['verbose']:
                    print("\tdoesn't exist. Skipping")
                continue
            df = pd.read_csv(results_file, sep='\t')

            # weight str for non-main weight types
            if weight_type != "":
                if main_weight_type and weight_type == main_weight_type[alg]:
                    weight_str = ""
                else:
                    weight_str = "%s-"%weight_type 
            if main_alg_opts and main_weight_type and alg_names and \
               alg in main_alg_opts and main_alg_opts[alg] in results_file and \
               (weight_type == "" or main_weight_type[alg] == weight_type):
                df['Algorithm'] = alg_names[alg]
            # include the alg-specific params in the name
            elif alg == "birgrank":
                df['Algorithm'] = 'birgrank %sa%s-t%s-m%s-l%s' % (weight_str, alpha, theta, mu, br_lambda)
            elif alg == "sinksource":
                df['Algorithm'] = 'sinksource %sa%s-maxi%s' % (weight_str, alpha, maxi)
            elif alg == "genemania":
                df['Algorithm'] = "genemania %stol%s" % (weight_str, tol)
            else:
                df['Algorithm'] = alg 

            if kwargs['verbose']:
                print("%d sp-goterm pairs" % (len(df['Algorithm'])))

            df_all = pd.concat((df_all, df))

    for alg in kwargs['algorithm']:
        alg_results = False
        for alg2 in df_all['Algorithm'].unique():
            if alg in alg2 or (alg_names and alg2 == alg_names[alg]):
                alg_results = True
                break
        if alg_results is False:
            "WARNING: no results found for algorithm %s" % (alg)
            time.sleep(2)

    return df_all


def results_overview(ev_code_results, measures=['fmax'], **kwargs):
    if 'main_alg_opts' in kwargs:
        print("main_alg_opts:")
        for alg in kwargs['algorithm']:
            print("\t%s: %s" % (alg, kwargs['main_alg_opts'][alg]))
    for version, exp_name in sorted(ev_code_results, key=lambda x: (x[1], x[0])):
        df_curr = ev_code_results[(version, exp_name)]
        print(version, exp_name)
        # limit the goterms to those that are also present for SinkSource(?)
        for measure in measures:
            for alg in sorted(df_curr['Algorithm'].unique()):
                print("\t%s: %0.3f \t\t(%d sp-goterm pairs)" % (alg, df_curr[df_curr['Algorithm'] == alg][measure].median(), len(df_curr[df_curr['Algorithm'] == alg][measure])))


def setup_optparse():
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    # general parameters
    group = OptionGroup(parser, 'Main Options (all can be comma-delimited)')
    group.add_option('','--version',type='string', default="2018_06-seq-sim-e0_1",
                     help="Version of the PPI to run. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    #group.add_option('-N','--net-file',type='string',
    #                 help="Network file to use. Can specify one per version. Default is the version's default network")
    group.add_option('-A', '--algorithm', default='sinksource',
                     help="Algorithm for which to get predictions. Default is 'sinksource' of them. Options: '%s'" % ("', '".join(alg_utils.ALGORITHMS)))
    #group.add_option('', '--exp-name', type='string', default='expc-50-1000-bp-use-neg',
    #                 help="Experiment name to used when running algorithm")
    # the experiment name is built with the following:
    # <ev_codes>-<goterm_size_range><eval_ev_codes>-<string_nets>-<h>-<use_neg>-<keep_ann>
    # example1: expc-50-1000-bp-use-neg
    # example2: expc-comp-50-1000iea-core-bp-use-neg-keep-ann
    group.add_option('', '--ev-codes', default='expc',
                     help="Evidence codes which define the pos/neg examples. Default='expc'")
    group.add_option('', '--size-range', default='50-1000',
                     help="GO term size range. Default='50-1000'")
    group.add_option('', '--eval-ev-codes', default='',
                     help="Evidence codes used to evaluate (if different than <ev-codes>. Default=''")
    group.add_option('', '--string-nets', default='core',
                     help="String networks used ('all', 'nontransferred', or 'core'). This will be left blank if string isn't in the version. Default='core'")
    group.add_option('-H', '--hierarchy', default='bp',
                     help="GO hierarchy. Either 'bp' or 'mf' or 'bp,mf'. Default='bp'")
    group.add_option('', '--pos-neg-str', default='-use-neg',
                     help="Either '', '-use-neg', '-non-pos-neg', '-use-neg-keep-ann', or '-non-pos-neg-keep-ann'. Default='-use-neg'")
    parser.add_option_group(group)

    add_alg_opts(parser)
    #run_algs.add_string_opts(parser)

    return parser


def add_alg_opts(parser):
    # parameters for running algorithms
    group = OptionGroup(parser, 'Algorithm options (all can be comma-delimited)')
    # should've put the string-nets and weight-type together...
    group.add_option('', '--weight-type', default="swsn",
                     help="Method used to integrate networks, or setup weights. Can be either '', 'unw', 'gm2008', or 'swsn'. Default='swsn' (for versions that include STRING)")
    group.add_option('-a', '--alpha', default="0.95,0.9",
                     help="Alpha insulation parameter, or restart parameter for BirgRank. Default=0.95")
    group.add_option('', '--eps', default="0.0,0.0001",
                     help="Stopping criteria for SinkSource or BirgRank. Default='0.0,0.0001'")
    group.add_option('', '--max-iters', default="10,20,1000",
                     help="Maximum # of iterations for SinkSource or BirgRank. Default=10,20,1000")
    group.add_option('', '--tol', default="1e-5",
                     help="Tolerance for convergance for the Scipy Sparse Conjugate Gradient solver (for GeneMANIA). Default=1e-05")
    group.add_option('', '--theta', default="0.5",
                     help="BirgRank: (1-theta) percent of Rtrain used in seeding vectors. Default=0.5")
    group.add_option('', '--mu', default="0.5",
                     help="BirgRank: (1-mu) percent of random walkers diffuse from G via Rtrain to H. Default=0.5")
    group.add_option('', '--br-lambda', default="0.01",
                     help="BirgRank: (1-lambda) percent of random walkers which diffuse downwards within H. Default=0.01")
    parser.add_option_group(group)

    # TODO 
#    # aptrank options
#    group = OptionGroup(parser, 'AptRank options')
#    group.add_option('', '--apt-k',
#                     help="Markov Chain iterations. Default=8")
#    group.add_option('', '--apt-s', 
#                     help="Number of shuffles. Default=5")
#    group.add_option('', '--apt-t', 
#                     help="Split percentage. Default=.5")
#    group.add_option('', '--diff-type', default="twoway",
#                     help="Diffusion type: oneway (G to H) or twoway (G to H and H to G). Default=twoway")
#    parser.add_option_group(group)


def validate_opts(opts):
    kwargs = vars(opts)

    #try:
    for opt in ['version', 'algorithm', 'ev_codes', 'eval_ev_codes', 
            'string_nets', 'hierarchy', 'weight_type']:
        kwargs[opt] = kwargs[opt].split(',')
    for opt in ['alpha', 'eps', 'tol', 'theta', 'mu', 'br_lambda']:
        kwargs[opt] = [float(param) for param in kwargs[opt].split(',')]
    for opt in ['max_iters']:
        kwargs[opt] = [int(param) for param in kwargs[opt].split(',')]
    #except:
    #    print("")

    for version in kwargs['version']:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)

    #if kwargs.algorithm not in f_settings.ALGORITHM_OPTIONS:
    #    print "--algorithm %s is not a valid algorithm name" % (kwargs.algorithm)
    #    sys.exit(1)

    # TODO
    for alg in kwargs['algorithm']:
        if alg not in alg_utils.ALGORITHMS:
            print("ERROR: '%s' not a valid algorithm name. Algorithms are: '%s'." % (alg, ', '.join(alg_utils.ALGORITHMS)))
            sys.exit(1)

    if kwargs['alg_plot_opts'] is not None:
        # set the gobal parameters
        global MAIN_WEIGHT_TYPE, MAIN_ALG_OPTS, ALG_NAMES
        alg_opts = json.load(open(kwargs['alg_plot_opts']))
        MAIN_WEIGHT_TYPE, MAIN_ALG_OPTS, ALG_NAMES = alg_opts['main_weight_type'], alg_opts['main_alg_opts'], alg_opts['alg_names']

    return kwargs


def add_plot_opts(parser):
    group = OptionGroup(parser, 'Plotting options')
    group.add_option('', '--forceplot', action="store_true", default=False,
                     help="Option to overwrite previous plots")
    group.add_option('', '--for-pub', action="store_true", default=False,
                     help="Option to remove headers and overall make plot nicer for publication")
    group.add_option('', '--postfix', default='',
                     help="str to add to end of filename")
    group.add_option('', '--alg-plot-opts', 
                     help="JSON file containing the parameters to use for each algorithm ('main_alg_opts'), weighting type ('main_weight_type'), and name ('alg_names').")
    parser.add_option_group(group)

    group = OptionGroup(parser, 'Additional options')
    group.add_option('', '--unweighted', action="store_true", default=False,
                     help="Option to ignore edge weights when running algorithms. Default=False (weighted)")
    group.add_option('-l', '--sinksourceplus-lambda', type=float, 
                     help="lambda parameter to specify the weight connecting the unknowns to the negative 'ground' node. Default=None")
    group.add_option('', '--verbose', action="store_true", default=False,
                     help="Print additional info about running times and such")
    parser.add_option_group(group)

    # TODO make this an option
    #group.add_option('', '--keep-ann', action='store_true', default=False,
    #                 help="Don't leave out annotations when running the algorithms" +
    #                 "TODO allow for the option to run and evaluate everything together")


def parse_args(args):
    parser = setup_optparse()
    add_plot_opts(parser)

    (opts, args) = parser.parse_args(args)
    kwargs = validate_opts(opts) 

    return kwargs


if __name__ == "__main__":
    kwargs = parse_args(sys.argv)
    versions = kwargs.pop('version')

    main(versions, **kwargs)
