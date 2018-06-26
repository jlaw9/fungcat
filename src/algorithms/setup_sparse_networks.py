
# Script to setup network and annotations files as sparse matrices
# also used for weighting the networks with the
# Simultaneous Weight with Specific Negatives (SWSN) method

import os, sys
from optparse import OptionParser, OptionGroup
from collections import defaultdict
#import matlab.engine
#from optparse import OptionParser
sys.path.append('src')
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils.file_utils as utils
import fungcat_settings as f_settings
from string_split_to_species import full_column_names, \
    STRING_NETWORKS, NON_TRANSFERRED_STRING_NETWORKS, \
    CORE_STRING_NETWORKS
import algorithms.gain_scipy.alg_utils as alg_utils
#from findKernelWeights import findKernelWeights
from algorithms.weight_networks.combineNetworksSWSN import combineNetworksSWSN
import networkx as nx
#import pandas as pd
import numpy as np
from scipy.io import savemat, loadmat
from scipy import sparse
from tqdm import tqdm
#import igacat.go_term_prediction_examples.go_term_prediction_examples as go_examples


def parse_args(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    # general parameters
    group = OptionGroup(parser, 'Main Options')
    group.add_option('','--version',type='string', action='append',
            help="Version of the PPI to run. Can specify multiple versions and they will run one after the other. " +
            "Default = '2017_10-string'" +
            "Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    group.add_option('', '--pos-neg-file', type='string', action='append',
            help="File containing positive and negative examples for each GO term. Required")
    group.add_option('-o', '--out-pref-net', type='string',
            help="Output prefix for the network file to create. " +
                     "Default: %s" % ('TODO'))
    group.add_option('', '--out-pref-ann', type='string',
            help="Output prefix for the annotations file to create. " +
                     "Default: %s" % ('TODO'))
    parser.add_option_group(group)

    # options to limit the # of GO terms or species
    group = OptionGroup(parser, 'Limitting options')
    group.add_option('', '--only-functions', type='string',
            help="Run GAIN using only the functions in a specified " +
            "file. If not specified, all functions will be used.")
    group.add_option('-G', '--goterm', type='string', action="append",
            help="Specify the GO terms to use (should be in GO:00XX format)")
    group.add_option('-T', '--taxon', type='string', action="append",
            help="Specify the species taxonomy ID to use.")
    parser.add_option_group(group)

    # parameters for STRING networks
    group = OptionGroup(parser, 'STRING options')
    group.add_option('', '--only-combined', action="store_true", default=False,
            help="Use only the STRING combined network: \n\tcombined_score")
    group.add_option('', '--only-core', action="store_true", default=False,
            help="Use only the 6 core networks: \n\t%s" % (', '.join(CORE_STRING_NETWORKS)))
    group.add_option('', '--non-transferred', action="store_true", default=False,
            help="Use all non-transferred networks: \n\t%s" % (', '.join(NON_TRANSFERRED_STRING_NETWORKS)))
    group.add_option('', '--all-string', action="store_true", default=False,
            help="Use all individual 13 STRING networks: \n\t%s" % (', '.join(STRING_NETWORKS)))
    group.add_option('-S', '--string-networks', type='string', 
            help="Comma-separated list of string networks to use. " +
                 "If specified, other STRING options will be ignored.")
    parser.add_option_group(group)

    # additional parameters
    group = OptionGroup(parser, 'Additional options')
    group.add_option('', '--weight-SWSN', type='string', 
            help="Weight and combine the networks using the" +
                 "Simultaneous Weights with Specific Negatives method." +
                 "Must specify a prefix for the network")
    group.add_option('', '--forcenet', action="store_true", default=False,
                     help="Force re-building network matrix from scratch")
    #group.add_option('', '--verbose', action="store_true", default=False,
    #                 help="Print additional info about running times and such")
    parser.add_option_group(group)

    (opts, args) = parser.parse_args(args)

    if opts.pos_neg_file is None:
        print("--pos-neg-file required")
        sys.exit(1)

    for f in opts.pos_neg_file:
        if not os.path.isfile(f):
            print("ERROR: --pos-neg-file '%s' not found" % (f))
            sys.exit(1)

    if opts.version is None:
        opts.version = ["2017_10-string"]

    for version in opts.version:
        if version not in f_settings.ALLOWEDVERSIONS:
            print("ERROR: '%s' not an allowed version. Options are: %s." % (version, ', '.join(f_settings.ALLOWEDVERSIONS)))
            sys.exit(1)

    return opts


def main(versions, pos_neg_files, goterms=None, taxons=None,
         out_pref_net=None, out_pref_ann=None,
         string_networks=["combined_score"],
         string_cutoff=f_settings.STRING_CUTOFF,
         weight_swsn=None, forcenet=False,
):

    for version in versions:
        INPUTSPREFIX, RESULTSPREFIX, network_file, selected_strains = f_settings.set_version(version)

        if out_pref_net is None:
            out_pref_net = '%s/sparse_nets/' % (INPUTSPREFIX)
        utils.checkDir(os.path.dirname(out_pref_net))
        if out_pref_ann is None:
            out_pref_ann = '%s/sparse_nets/' % (INPUTSPREFIX)
        utils.checkDir(os.path.dirname(out_pref_ann))

        if taxons is not None:
            for taxon in taxons:
                out_pref = "%s%s-" % (out_pref_net, taxon)
                sparse_networks, network_names, nodes = create_sparse_net_file(
                    version, out_pref, taxon=taxon,
                    string_nets=string_networks, string_cutoff=string_cutoff,
                    forcenet=forcenet)

                out_pref = "%s%s-" % (out_pref_ann, taxon)
                ann_matrix, goids = create_sparse_ann_file(
                    out_pref, pos_neg_files, goterms, nodes,
                    selected_strains=selected_strains, taxon=taxon)

                if weight_swsn:
                    out_file = "%s%s-%d-nets-combined-SWSN.npz" % (
                        weight_swsn, taxon, len(sparse_networks))
                    weight_SWSN(ann_matrix, sparse_networks,
                                net_names=network_names, out_file=out_file,
                                nodes=nodes)
        else:
            sparse_networks, network_names, nodes = create_sparse_net_file(
                version, out_pref_net, selected_strains=selected_strains,
                string_nets=string_networks, string_cutoff=string_cutoff,
                forcenet=forcenet)

            # now write the annotations
            ann_matrix, goids = create_sparse_ann_file(
                out_pref_ann, pos_neg_files, goterms, nodes,
                selected_strains=selected_strains)

            if weight_swsn:
                out_file = "%s%d-nets-combined-SWSN.npz" % (
                    weight_swsn, len(sparse_networks))
                weight_SWSN(ann_matrix, sparse_networks,
                            net_names=network_names, out_file=out_file,
                            nodes=nodes)


def create_sparse_net_file(
        version, out_pref, selected_strains=None, taxon=None, string_nets=None,
        string_cutoff=f_settings.STRING_CUTOFF, forcenet=False):

    net_files = []
    string_net_files = []
    if taxon is not None:
        network_file = f_settings.STRING_TAXON_UNIPROT_FULL % (taxon, taxon, string_cutoff)
        #ann_file = f_settings.FUN_FILE % (self.taxon, self.taxon)
        string_net_files.append(network_file)

    else:
        # use all of the networks available to this version by default
        if 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[version]:
            # load the seq sim network for this version
            net_files.append(f_settings.SEQ_SIM_NETWORKS[version])
        if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[version]:
            # and all of the string networks available
            for taxon in selected_strains:
                net_file = f_settings.STRING_TAXON_UNIPROT_FULL % (taxon, taxon, string_cutoff)
                string_net_files.append(net_file)

    # the node IDs should be the same for each of the networks,
    # so no need to include the # in the ids file
    num_networks = len(net_files) + len(string_nets)
    sparse_nets_file = "%s%d-sparse-nets.mat" % (out_pref, num_networks)
    node_ids_file = "%snode-ids.txt" % (out_pref)
    net_names_file = "%s%d-net-names.txt" % (out_pref, num_networks)
    if forcenet is False \
       and os.path.isfile(sparse_nets_file) and os.path.isfile(node_ids_file) \
       and os.path.isfile(net_names_file):
        # read the files
        print("\treading sparse nets from %s" % (sparse_nets_file))
        sparse_networks = loadmat(sparse_nets_file)['Networks'][0]
        print("\treading node ids file from %s" % (node_ids_file))
        nodes = utils.readItemList(node_ids_file, 1)
        print("\treading network_names from %s" % (net_names_file))
        network_names = utils.readItemList(net_names_file, 1)

    else:
        print("\tcreating sparse nets and writing to %s" % (sparse_nets_file))
        sparse_networks, network_names, nodes = setup_sparse_networks(
            net_files=net_files, string_net_files=string_net_files, string_nets=string_nets)

        # now write them to a file
        write_sparse_net_file(
            sparse_networks, sparse_nets_file, network_names,
            net_names_file, nodes, node_ids_file)

    return sparse_networks, network_names, nodes


def write_sparse_net_file(
        sparse_networks, out_file, network_names,
        net_names_file, nodes, node_ids_file):
    #out_file = "%s/%s-net.mat" % (out_dir, version)
    # save this graph into its own matlab file because it doesn't change by GO category
    print("\twriting sparse networks to %s" % (out_file))
    mat_networks = np.zeros(len(sparse_networks), dtype=np.object)
    for i, net in enumerate(network_names):
        mat_networks[i] = sparse_networks[i]
    savemat(out_file, {"Networks":mat_networks}, do_compression=True)

    print("\twriting node2idx labels to %s" % (node_ids_file))
    with open(node_ids_file, 'w') as out:
        out.write(''.join(["%s\t%s\n" % (n, i) for i, n in enumerate(nodes)]))

    print("\twriting network names, which can be used to figure out " +
          "which network is at which index in the sparse-nets file, to %s" % (net_names_file))
    with open(net_names_file, 'w') as out:
        out.write(''.join(["%s\n" % (n) for n in network_names]))

    # TODO can't save multiple networks to a npz file.
    # for now, just write the .mat file
    #out_file = out_file.replace('.mat', '.npz')
    #print("\twriting sparse network to %s" % (out_file))
    #sparse.save_npz(out_file, sparse_matrix)

    #if os.path.isfile(prots_file) and os.path.isfile(out_file):
    #    print("Skipping mapping nodes to ids. %s exists" % (out_file))
    #    print("Reading the prot ordering/node2idx mapping from %s" % (prots_file))
    #    prots = utils.readItemList(prots_file, 1)
    #    node2idx = {prot: int(idx) for prot, idx in utils.readColumns(prots_file, 1, 2)}
    #    normalized_nets = loadmat(out_file)['Networks'][0] 
    #else:
    #    network_file = "inputs/%s/%s-net.txt" % (version, version)


def setup_sparse_networks(net_files=[], string_net_files=[], string_nets=[]):
    """
    Function to setup networks as sparse matrices 
    *net_files*: list of networks for which to make into a sparse
        matrix. The name of the file will be the name of the sparse matrix
    *string_net_files*: List of string files containing all 14 STRING network columns
    *string_nets*: List of STRING network column names for which to make a sparse matrix. 

    *returns*: List of sparse networks, list of network names, 
        list of proteins in the order they're in in the sparse networks
    """

    network_names = []
    # TODO build the sparse matrices without using networkx
    # I would need to make IDs for all proteins to ensure the IDs and
    # dimensions are the same for all of the matrices
    G = nx.Graph()
    for net_file in tqdm(net_files):
        name = os.path.basename(net_file)
        network_names.append(name)
        tqdm.write("Reading network from %s. Giving the name %s" % (net_file, name))
        with open(net_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue
                #u,v,w = line.rstrip().split('\t')[:3]
                line = line.rstrip().split('\t')
                u,v,w = line[:3]
                G.add_edge(u,v,**{name:float(w)})

    network_names += string_nets
    print("Reading %d STRING networks" % len(string_net_files))
    # for now, group all of the species string networks into a
    # massive network for each of the string_nets specified 
    for string_net_file in tqdm(string_net_files):
        tqdm.write("Reading network from %s" % (string_net_file))
        with open(string_net_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue
                #u,v,w = line.rstrip().split('\t')[:3]
                line = line.rstrip().split('\t')
                u,v = line[:2]
                attr_dict = {}
                for net in string_nets:
                    attr_dict[net] = float(line[full_column_names[net]-1])
                # if the edge already exists, it will be added with
                # the new attributes
                G.add_edge(u,v,**attr_dict)
    print("\t%d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))

    print("\trelabeling node IDs with integers")
    G, node2idx, idx2node = alg_utils.convert_nodes_to_int(G)
    # keep track of the ordering for later
    nodes = [idx2node[n] for n in sorted(idx2node)]

    print("\tconverting graph to sparse matrices")
    #sparse_networks = np.zeros(len(net_files)+len(string_net_files), dtype=np.object)
    sparse_networks = []
    for i, net in enumerate(tqdm(network_names)):
        # default nodelist order is G.nodes()
        sparse_matrix = nx.to_scipy_sparse_matrix(G, nodelist=sorted(idx2node), weight=net)
        # convert to float otherwise matlab won't parse it correctly
        # see here: https://github.com/scipy/scipy/issues/5028
        sparse_matrix = sparse_matrix.astype(float) 
        sparse_networks.append(sparse_matrix)

    # save this graph into its own matlab file because it doesn't change by GO category
    #print("\twriting sparse networks to %s" % (out_file))
    #savemat(out_file, {"Networks":sparse_networks})

    # TODO can't save multiple networks to a npz file.
    # for now, just write the .mat file
    #out_file = out_file.replace('.mat', '.npz')
    #print("\twriting sparse network to %s" % (out_file))
    #sparse.save_npz(out_file, sparse_matrix)

    return sparse_networks, network_names, nodes


def create_sparse_ann_file(out_pref, pos_neg_files, goids, prots,
                           selected_strains=None, taxon=None):
    # setup the annotations
    ann_matrix, goids = setup_sparse_annotations(
        pos_neg_files, goids, prots, selected_species=selected_strains, taxon=taxon)

    out_file = "%s%d-annotations.mat" % (out_pref, len(goids))
    # convert the GO DAG to a sparse matrix, while maintaining the order of goids so it matches with the annotation matrix
    #dag_matrix = nx.to_scipy_sparse_matrix(G, nodelist=goids, weight=None)
    #print("\twriting sparse annotations to %s" % (out_file))
    #savemat(out_file, {'annotations': ann_matrix})
    # also save it as a scipy file
    out_file = out_file.replace('.mat', '.npz')
    print("\twriting sparse annotations to %s" % (out_file))
    sparse.save_npz(out_file, ann_matrix)

    # write this order of the GO IDs to a file so they can easily be accessed (to get the rows of the sparse matrix)
    # without having to perform these operations again
    goids_file = "%s%d-goids.txt" % (out_pref, len(goids))
    print("Writing %s" % (goids_file))
    with open(goids_file, 'w') as out:
        out.write(''.join("%s\t%d\n" % (goid, i) for i, goid in enumerate(goids)))
    #goid2idx = {goid: i for i, goid in enumerate(goids)}

    return ann_matrix, goids


def setup_sparse_annotations(pos_neg_files, goterms, prots,
                             selected_species=None, taxon=None):
    """
    
    *returns*: 1) A matrix with goterm rows, protein/node columns, and
        1,0,-1 for pos,unk,neg values
        2) List of goterms in the order in which they appear in the matrix
    """
    # TODO figure out a better way to keep track of the file
    #num_goterms = 6
    #h = "bp"
    #ann_file = "%s/%s-%d-annotations.npz" % (out_dir, h, num_goterms)
    #
    #if os.path.isfile(ann_file):
    #    print("\nReading annotations from %s" % (ann_file))
    #    ann_matrix = sparse.load_npz(ann_file)
    #    print("\t%d goterms, %d pos/neg annotations" %
    #          (ann_matrix.shape[0], len(ann_matrix.data)))
    #    goids_file = "%s/%s-%d-goids.txt" % (out_dir, h, num_goterms)
    #    goids = utils.readItemList(goids_file, 1)
    #else:
    print("\nSetting up annotation matrix")
    #selected_species_file = f_settings.VERSION_SELECTED_STRAINS[version]
    #selected_species = utils.readDict(selected_species_file, 1, 2)

    # this file just has the direct GO annotations (not propagated)
    #gaf_file = "/data/inputs/goa/taxon/19-strains-goa.gaf"
    #ev_codes = "rem-neg-iea"
    #h = "bp"
    ## TODO make a pos_neg_file for each species
    #pos_neg_files = [
    #        "inputs/pos-neg/%s/pos-neg-%s-50-list.tsv" % (ev_codes, h),] 
    ##        "inputs/pos-neg/%s/pos-neg-mf-50-list.tsv" % (ev_codes),]
    ## TODO build the matrix while parsing the file
    goid_pos, goid_neg = alg_utils.parse_pos_neg_files(pos_neg_files, goterms=goterms) 
    node2idx = {prot: i for i, prot in enumerate(prots)}

    # limit it to the current taxon
    if taxon is not None:
        print("Getting species of each prot from %s" % (f_settings.UNIPROT_TO_SPECIES))
        #print("Limiting the prots to those for taxon %s (%s)" % (taxon, selected_species[taxon]))
        print("Limiting the prots to those for taxon %s" % (taxon))
        # for each of the 19 species, leave out their annotations 
        # and see how well we can retrieve them 
        uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1,2)
        # also build the reverse
        species_to_uniprot = defaultdict(set)
        for p in uniprot_to_species:
            species_to_uniprot[uniprot_to_species[p]].add(p)
        if taxon not in species_to_uniprot:
            print("Error: taxon ID '%d' not found" % (taxon))
            sys.exit()
        taxon_prots = species_to_uniprot[taxon]
        # also limit the proteins to those in the network
        print("\t%d prots for taxon %s. Limiting to the %d in the network" % (len(taxon_prots), taxon, len(node2idx)))
        taxon_prots = taxon_prots & set(node2idx.keys())
        goid_pos = {goid: prots & taxon_prots for goid, prots in goid_pos.items()}
        goid_neg = {goid: prots & taxon_prots for goid, prots in goid_neg.items()}
        goids_to_remove = set()
        for goid in goid_pos:
            if len(goid_pos[goid]) == 0:
                goids_to_remove.add(goid)
        if len(goids_to_remove) > 0:
            print("Warning: %d goterms have 0 annotations. Removing them." % (len(goids_to_remove)))
            for goid in goids_to_remove:
                goid_pos.pop(goid)
                goid_neg.pop(goid)

    goids = sorted(goid_pos.keys())

#def build_sparse_matrix(data, rows, cols):
    # rather than explicity building the matrix, use the indices to build a coordinate matrix
    # rows are prots, cols are goids
    i_list = []
    j_list = []
    data = []
    num_pos = 0
    num_neg = 0
    # limit the annotations to the proteins which are in the networks
    prots_set = set(prots)
    for j, goid in enumerate(goids):
        for prot in goid_pos[goid] & prots_set:
            i_list.append(node2idx[prot])
            j_list.append(j)
            data.append(1)
            num_pos += 1
        for prot in goid_neg[goid] & prots_set:
            i_list.append(node2idx[prot])
            j_list.append(j)
            data.append(-1)
            num_neg += 1
    print("\t%d annotations. %d positive, %d negatives" % (len(data), num_pos, num_neg))

    # convert it to a sparse matrix 
    print("Building a sparse matrix of annotations")
    ann_matrix = sparse.coo_matrix((data, (i_list, j_list)), shape=(len(prots), len(goids)), dtype=float).tocsr()
    print("\t%d pos/neg annotations" % (len(ann_matrix.data)))
    ann_matrix = ann_matrix.transpose()

    return ann_matrix, goids

#    # convert the GO DAG to a sparse matrix, while maintaining the order of goids so it matches with the annotation matrix
#    #dag_matrix = nx.to_scipy_sparse_matrix(G, nodelist=goids, weight=None)
#    out_file = "%s/%s-%d-annotations.mat" % (out_dir, h, len(goids))
#    print("\twriting %s" % (out_file))
#    savemat(out_file, {'annotations': ann_matrix})
#    # also save it as a scipy file
#    out_file = out_file.replace('.mat', '.npz')
#    print("\twriting sparse annotations to %s" % (out_file))
#    sparse.save_npz(out_file, ann_matrix)
#    #savemat(out_file, {'R':annotation_matrix, 'H': dag_matrix})
#
#    ann_matrix = ann_matrix[0:20]
#


def weight_SWSN(ann_matrix, sparse_nets, net_names=None, out_file=None, nodes=None):
    """
    """

    # normalize the networks
    print("Normalizing the networks")
    normalized_nets = []
    for net in sparse_nets:
        normalized_nets.append(_net_normalize(net))
    print("Weighting networks for %d different GO terms" % (ann_matrix.shape[0]))
    print("Running simultaneous weights with specific negatives")
    # try running it now
    #for i in range(ann_matrix.shape[0]):
    #    y = ann_matrix[i].toarray()[0]
    #    print("\tgoid %s: %d positives, %d negatives" % (goids[i], len(np.where(y > 0)[0]), len(np.where(y < 0)[0])))
    #    #alpha, indices = findKernelWeights(y, normalized_nets)
    #    #print("\tString networks: %s\n" % (', '.join(STRING_NETWORKS[x] for
    #    #                                           x in indices)))
    alpha, indices = combineNetworksSWSN(ann_matrix, normalized_nets) 
    if net_names is not None:
        print("\tnetworks chosen: %s" % (', '.join([net_names[i] for i in indices])))

    # now add the networks together with the alpha weight applied
    combined_network = alpha[0]*sparse_nets[indices[0]]
    for i in range(1,len(alpha)):
        combined_network += alpha[i]*sparse_nets[indices[i]] 

    if out_file is not None:
        utils.checkDir(os.path.dirname(out_file))
        print("\twriting combined network to %s" % (out_file))
        sparse.save_npz(out_file, combined_network)
        # also write the node ids so it's easier to access
        # TODO figure out a better way to store this
        node2idx_file = out_file + "-node-ids.txt"
        print("\twriting node ids to %s" % (node2idx_file)) 
        with open(node2idx_file, 'w') as out:
            out.write(''.join("%s\t%s\n" % (n, i) for i, n in enumerate(nodes)))

        # write the alpha/weight of the networks as well
        net_weight_file = out_file + "-net-weights.txt"
        print("\twriting network weights to %s" % (net_weight_file)) 
        with open(net_weight_file, 'w') as out:
            out.write(''.join("%s\t%s\n" % (net_names[idx], str(alpha[i])) for i, idx in enumerate(indices)))

    return combined_network


# this was mostly copied from the deepNF preprocessing script
def _net_normalize(X):
    """ 
    Normalizing networks according to node degrees.
    """
    if X.min() < 0:
        print("### Negative entries in the matrix are not allowed!")
        X[X < 0] = 0 
        print("### Matrix converted to nonnegative matrix.")
    # for now assume the network is symmetric
    #if (X.T != X).all():
    #    pass
    #else:
    #    print("### Matrix not symmetric.")
    #    X = X + X.T - np.diag(np.diag(X))
    #    print("### Matrix converted to symmetric.")

    # normalizing the matrix
    deg = X.sum(axis=1).A.flatten()
    deg = np.divide(1., np.sqrt(deg))
    deg[np.isinf(deg)] = 0 
    # sparse matrix function to make a diagonal matrix
    D = sparse.spdiags(deg, 0, X.shape[0], X.shape[1], format="csr")
    X = D.dot(X.dot(D))

    return X


def run():
    #versions = ["2017_10-seq-sim", "2017_10-seq-sim-x5-string"]
    opts = parse_args(sys.argv)
    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 
    print("%d GO terms from only_functions_file and/or specified GO terms" % (0 if goterms is None else len(goterms)))

    #goid_pos, goid_neg = alg_utils.parse_pos_neg_files(opts.pos_neg_file) 

    # setup the selection of string networks 
    if opts.string_networks:
        string_networks = opts.string_networks.split(',')
        for net in string_networks:
            if net not in STRING_NETWORKS:
                print("ERROR: STRING network '%s' not one of the" +
                      "available choices which are: \n\t%s" % (net, ', '.join(STRING_NETWORKS)))
                sys.exit(1)
    elif opts.only_combined:
        string_networks = ['combined_score']
    elif opts.only_core:
        string_networks = CORE_STRING_NETWORKS
    elif opts.non_transferred:
        string_networks = NON_TRANSFERRED_STRING_NETWORKS
    elif opts.all_string:
        string_networks = STRING_NETWORKS

    main(opts.version, opts.pos_neg_file, goterms=goterms, taxons=opts.taxon,
         out_pref_net=opts.out_pref_net, out_pref_ann=opts.out_pref_ann,
         string_networks=string_networks, string_cutoff=f_settings.STRING_CUTOFF,
         weight_swsn=opts.weight_SWSN, forcenet=opts.forcenet
    )


if __name__ == "__main__":
    run()
