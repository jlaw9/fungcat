
from optparse import OptionParser, OptionGroup
from collections import defaultdict
import os
import sys 
from tqdm import tqdm
import itertools
sys.path.append('src')
import utils.file_utils as utils
import fungcat_settings as f_settings
#import algorithms.gain_scipy.run_algs as run_algs
import algorithms.gain_scipy.alg_utils as alg_utils
# import plotting libraries
import matplotlib
matplotlib.use('Agg') # To save files remotely.  Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np


# tried to be fancy :P
# colors: https://coolors.co/ef6e4a-0ec9aa-7c9299-5d88d3-96bd33
my_palette = ["#EF6E4A", "#0EC9AA", "#7C9299", "#5D88D3", "#96BD33", "#937860", "#efd2b8"]

# these are the main options used to determine which parameters to use for the algorithm
# (unless the plot is specifically comparing parameters, then all will be used)
# for now, just check if the entire string is in the results file. If it is, then use it
MAIN_ALG_OPTS = {
    'sinksource': 'a0_95-eps0_0-maxi10',
    'genemania': "a0_95-eps0_0-maxi20-tol1e-05",  # TODO remove a, eps and maxi from the output file name
    'birgrank': "a0_9-eps0_0001-maxi1000-t0_5-m0_5-l0_01",
    'localplus': "a0_95-eps0_0-maxi20",  # TODO remove a, eps and maxi from the output file name
    }

ALG_NAME = {
    'localplus': 'Local', 'local': 'Local-',
    'sinksource': 'SinkSource', 'sinksourceplus': 'SinkSource+',
    'genemania': 'GeneMANIA',  
    'birgrank': 'BirgRank', 'aptrank': 'AptRank OneWay',
    'sinksource-squeeze': 'SinkSource Squeeze', 
    }


def main(versions, **kwargs):
    #weight_str = 'swsn'
    #for version in kwargs['version']:
    #    get_results(version, weight_str, **kwargs)

    get_results(versions, **kwargs)


def get_results(versions, **kwargs):
    for version in versions:
        # for now, just use SWSN
        weight_str = "-swsn" if 'string' in version else ''

        for ev_codes, eval_ev_codes, string, h in itertools.product(
                *[kwargs['ev_codes'], kwargs['eval_ev_codes'], kwargs['string_nets']+[''], kwargs['hierarchy']]):
            exp_name = "%s-%s%s%s-%s%s" % (ev_codes, kwargs['size_range'], eval_ev_codes, '-%s'%string if string != '' else '', h, kwargs['pos_neg_str'])

            if 'string' not in version and \
               ('core' in exp_name or 'all' in exp_name or 'nontransferred' in exp_name):
                continue
            elif 'string' in version and \
                 ('core' not in exp_name and 'all' not in exp_name and 'nontransferred' not in exp_name):
                continue

            df_all = get_results_alg_params(version, exp_name, weight_str=weight_str, **kwargs)
            print(df_all.head())

            # TODO make an option for this
            cutoffs = [(10,1000)]
            if 's200' in version:
                cutoffs = [(10,5000)]
            #cutoffs_data = split_results_by_cutoffs(df_all, cutoffs, exp_goid_prots, comp_goid_prots, iea_goid_prots, overall_goid_prots=overall_goid_prots)
            # store this for later
            #ev_code_results[(version, ev_codes, eval_ev_codes, keep_ann, h)] = (cutoffs_data, out_dir)


#def get_results(version, weight_str, **kwargs):
def get_results_alg_params(version, exp_name, weight_str='', **kwargs):
    """ Get the fmax results for a given set of algorithm parameters
    *weight_str*: can be either '', '-unw', '-gm2008', or '-swsn'
    """
    df_all = pd.DataFrame()
    for alg in kwargs['algorithm']:
        for alpha, theta, mu, br_lambda, eps, maxi, tol in itertools.product(
                *[kwargs['alpha'], kwargs['theta'], kwargs['mu'], kwargs['br_lambda'], kwargs['eps'], kwargs['max_iters'], kwargs['tol']]):
            results_file = alg_utils.get_filepath(
                version, alg, exp_name, weight_str=weight_str, exp_type='loso', 
                alpha=alpha, theta=theta, mu=mu, br_lambda=br_lambda, eps=eps, max_iters=maxi, tol=tol,
                sinksourceplus_lambda=kwargs['sinksourceplus_lambda']) 

            if kwargs['verbose']:
                print("Reading results from %s" % (results_file))
            if not os.path.isfile(results_file):
                if kwargs['verbose']:
                    print("\tdoesn't exist. Skipping")
                continue
            df = pd.read_csv(results_file, sep='\t')

            if MAIN_ALG_OPTS[alg] in results_file:
                df['Algorithm'] = ALG_NAME[alg]
            else:
                df['Algorithm'] = ALG_NAME[alg]
            if kwargs['verbose']:
                print("%d sp-goterm pairs" % (len(df['Algorithm'])))

            df_all = pd.concat((df_all, df))

    return df_all


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
    for opt in ['version', 'algorithm', 'ev_codes', 'eval_ev_codes', 'string_nets', 'hierarchy']:
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

    return kwargs


def parse_args(args):
    parser = setup_optparse()

    group = OptionGroup(parser, 'Additional options')
    group.add_option('', '--unweighted', action="store_true", default=False,
                     help="Option to ignore edge weights when running algorithms. Default=False (weighted)")
    group.add_option('-l', '--sinksourceplus-lambda', type=float, 
                     help="lambda parameter to specify the weight connecting the unknowns to the negative 'ground' node. Default=None")
    group.add_option('', '--verbose', action="store_true", default=False,
                     help="Print additional info about running times and such")
    # TODO make this an option
    #group.add_option('', '--keep-ann', action='store_true', default=False,
    #                 help="Don't leave out annotations when running the algorithms" +
    #                 "TODO allow for the option to run and evaluate everything together")

    (opts, args) = parser.parse_args(args)
    kwargs = validate_opts(opts) 

    return kwargs


if __name__ == "__main__":
    kwargs = parse_args(sys.argv)
    versions = kwargs.pop('version')

    main(versions, **kwargs)

