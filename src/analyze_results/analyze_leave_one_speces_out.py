#! /usr/bin/python

# script to analyze the results of the leave-one-species-out evaluation
# specifically, compute the distance of left-out positives to positives,
#     left-out negatives to positives, left-out positives to negatives,
#     and left-out negatives to negatives


import sys, os
from collections import defaultdict
#from optparse import OptionParser
import networkx as nx
sys.path.append('src')
import utils.file_utils as utils
import fungcat_settings as f_settings
#import utils.graphspace.post_to_graphspace as gs
import algorithms.gain_scipy.alg_utils as alg_utils
import algorithms.gain_scipy.eval_leave_one_species_out as eval_sp
#from pandas import read_csv
from scipy import sparse
import numpy as np
from tqdm import tqdm
import fcntl


def main(version, exp_name, goid_pos, goid_neg,
         algorithms, net_file=None, taxon=None, unweighted=False,
         eval_pos=None, eval_neg=None):
    """
    *eval_pos/neg*: dictionaries with positive and negative examples to compute shortest paths to
    """

    if opts.non_pos_as_neg_eval is False:
        exp_name += "-use-neg" 
    else:
        exp_name += "-non-pos-neg" 

    # set the version for the UNIPROT_TO_SPECIES
    _, RESULTSPREFIX, NETWORK, _ = f_settings.set_version(version)

    use_scipy = False
    if opts.net_file is not None:
        use_scipy = True
        print("\tloading network from %s" % (opts.net_file))
        W, prots = alg_utils.setup_sparse_network(opts.net_file)
        node2idx = {n: i for i, n in enumerate(prots)}
        net_nodes = set(list(prots))
        # make sure to normalize it
        P = alg_utils.normalizeGraphEdgeWeights(W)
        # take the -log so smaller weights are better
        P.data = -np.log10(P.data)
    else:
        print("Reading network from %s" % (NETWORK))
        G = nx.Graph(name=version)
        #if unweighted is True:
        #    G = nx.read_edgelist(NETWORK, delimiter="\t", create_using=G)
        #else:
        G = nx.read_weighted_edgelist(NETWORK, delimiter="\t", create_using=G)
        print("\t%d nodes, %d edges" % (G.number_of_nodes(), G.number_of_edges()))
        net_nodes = set(G.nodes())

    out_dir = "outputs/viz/eval-species/distances/%s/%s" % (exp_name, version)
    utils.checkDir(out_dir)
    # distance from left-out positives to positives
    pos_pos_file = "%s/pos-pos%s.tsv" % (out_dir, "" if unweighted is False else "-unw")
    # distance from left-out positives to negatives
    pos_neg_file = "%s/pos-neg%s.tsv" % (out_dir, "" if unweighted is False else "-unw")

    print("Getting species of each prot from %s" % (f_settings.VERSION_UNIPROT_TO_SPECIES[version]))
    #print("Limiting the prots to those for taxon %s (%s)" % (taxon, selected_species[taxon]))
    print("Limiting the prots to those for taxon %s" % (taxon))
    # for each of the 19 species, leave out their annotations 
    # and see how well we can retrieve them 
    uniprot_to_species = utils.readDict(f_settings.VERSION_UNIPROT_TO_SPECIES[version], 1,2)
    # also build the reverse
    species_to_uniprot = defaultdict(set)
    for p in uniprot_to_species:
        species_to_uniprot[uniprot_to_species[p]].add(p)
    if taxon not in species_to_uniprot:
        print("Error: taxon ID '%d' not found" % (taxon))
        sys.exit()
    taxon_prots = species_to_uniprot[taxon]

    for goid in tqdm(goid_pos):
        # limit it to the current taxon
        #if taxon is not None:
        #print("\t%d prots for taxon %s." % (len(taxon_prots), taxon))
        non_taxon_pos_prots = set(goid_pos[goid]) - taxon_prots 
        non_taxon_neg_prots = set(goid_neg[goid]) - taxon_prots 
        if eval_pos is not None:
            if goid not in eval_pos or goid not in eval_neg:
                continue
            # if I'm evaluating COMP or ELEC, then compute the distance from those prots to the non-left-out prots
            pos_prots = set(eval_pos[goid]) & taxon_prots 
            neg_prots = set(eval_neg[goid]) & taxon_prots 
        else:
            pos_prots = set(goid_pos[goid]) & taxon_prots 
            neg_prots = set(goid_neg[goid]) & taxon_prots

        #print("\t%d annotated prots for %s" % (len(pos_prots), goid))
        # also limit the proteins to those in the network
        pos_prots = pos_prots & net_nodes
        neg_prots = neg_prots & net_nodes
        non_taxon_pos_prots = non_taxon_pos_prots & net_nodes
        non_taxon_neg_prots = non_taxon_neg_prots & net_nodes
        if len(pos_prots) < opts.num_test_cutoff or len(non_taxon_pos_prots) < opts.num_test_cutoff:
            print("\t%s: %d pos or %d non-taxon pos < %d num_test_cutoff. Skipping " % (goid, len(pos_prots), len(non_taxon_pos_prots), opts.num_test_cutoff))
            continue
        
        print("\t%s: %d pos, %d neg, %d non-taxon pos, %d non-taxon neg prots are in the network" % (
            goid, len(pos_prots), len(neg_prots), len(non_taxon_pos_prots), len(non_taxon_neg_prots)))

        # TODO for some reason, the shortest paths routine is giving incorrect results, and it's giving different results each time I run it
#        # now load the network as a sparse matrix
#        P, node2idx, idx2node = load_network(net_file=net_file, version_net_file=NETWORK, ss_lambda=None)
#        net_nodes = set(node2idx.keys())
#
#        # limit to nodes in the network
#        pos_prots = [node2idx[p] for p in pos_prots & net_nodes]
#        neg_prots = [node2idx[p] for p in neg_prots & net_nodes]
#        indices = pos_prots + neg_prots
#        if taxon is not None:
#            non_taxon_pos_prots = [node2idx[p] for p in non_taxon_pos_prots & net_nodes]
#            non_taxon_neg_prots = [node2idx[p] for p in non_taxon_neg_prots & net_nodes]
#            indices = pos_prots + non_taxon_pos_prots
#
#        # little sanity check:
#        #p1 = "P58619"
#        #p2 = "Q9HZ46"
#        #print("%s - %s: %s" % (p1, p2, P[node2idx[p1],node2idx[p2]]))
#        #sys.exit()
#
#        # find the shortest paths from the left-out positives to the positives of the other species
#        print("Finding the shortest paths for %d taxon prots and %d non-taxon prots (%d total)"
#              % (len(pos_prots), len(non_taxon_pos_prots), len(indices)))
#        shortest_paths = sparse.csgraph.shortest_path(P, directed=False, unweighted=unweighted, indices=indices) 
#        print(shortest_paths)
#        #for i, n in enumerate([pos_prots[0]]):
#        #    for j, n2 in enumerate(non_taxon_pos_prots):
#        #        print("shortest path from %s to %s: %s" % (n, n2, str(shortest_paths[i,j])))
#
#        # now write to a file / plot a histogram
#        out_file = "outputs/viz/eval-species/distances/%s-%s-%s%s.tsv" % (
#            version, goid, taxon, "" if unweighted is False else "-unw")
#        print("\tWriting to %s" % (out_file))
#        utils.checkDir(os.path.dirname(out_file))
#        with open(out_file, 'w') as out:
#            out.write(''.join("%s\t%s\t%s\n" % (
#                idx2node[n], idx2node[n2], shortest_paths[i,j])
#                for i, n in enumerate(pos_prots)
#                for j, n2 in enumerate(non_taxon_pos_prots)))

        # find the connected components
#        # first remove the edges with positives and negatives from the graph
#        # find the connected components,
#        # then add the pos/neg back in
#        #removed_edges = set()
#        #for u in non_taxon_neg_prots | non_taxon_pos_prots:
#        #    removed_edges.update(set([(u,v) for v in nx.all_neighbors(G, u)]))
#        G.remove_nodes_from(non_taxon_neg_prots | non_taxon_pos_prots)
#
#        components = nx.connected_components(G)
#
#        # now find the connected component around some of the positives
#        for i, cc in enumerate(sorted(components, key=len, reverse=True)):
#            pos_in_cc = set(cc) & pos_prots
#            if len(pos_in_cc) > 0 or len(cc) > 100:
#                print("%d positives (%s) in %d component (%d nodes)" % (len(pos_in_cc), ', '.join(pos_in_cc), i, len(cc)))

        goid_out_file = None
        # store the path lengths for each goid
        #goid_path_lengths = {}
        # TODO write a summary of the path lengths for each goid?
#        if len(goid_pos) <= 2:
#            start_end_sets = [
#                # left-out positives to themselves
#                (pos_prots, pos_prots, 'pos'),  
#                # find the shortest paths from the left-out positives to the positives of the other species
#                (pos_prots, non_taxon_pos_prots, 'pos-pos'),
#                # from the left-out positives to the negatives
#                (pos_prots, non_taxon_neg_prots, 'pos-neg'),
#                # from the left-out negatives to the negatives
#                (neg_prots, non_taxon_pos_prots, 'neg-pos')
#                ]
#
#            for start_prots, end_prots, ex_str in start_end_sets:
#                goid_out_file = "outputs/viz/eval-species/distances/%s-%s-%s%s-%s-nx.tsv" % (
#                    version, goid, taxon, "" if unweighted is False else "-unw", ex_str)
#                path_lengths = find_shortest_paths(G, start_prots, end_prots,
#                                                unweighted=True, out_file=goid_out_file)
#        else:
        # TODO set an option for this
        start_end_sets = [
            # find the shortest paths from the left-out positives to the positives of the other species
            (pos_prots, non_taxon_pos_prots, pos_pos_file),
            # from the left-out positives to the negatives
            (pos_prots, non_taxon_neg_prots, pos_neg_file),
        ]

        if use_scipy:
            #path_lengths = find_shortest_paths(W, start_prots, end_prots, unweighted=unweighted)
            # returns the distance of each source (row) to every other node in the graph
            start_prots = list(pos_prots)
            start_idx = [node2idx[p] for p in start_prots]
            dist_matrix = sparse.csgraph.shortest_path(P, directed=True, indices=start_idx)
            for _, end_prots, out_file in start_end_sets:
                end_prots = list(end_prots)
                if len(goid_pos) <= 2:
                    out_file = out_file.replace(".txt","%s-%s.txt" % (taxon, goid))
                end_idx = [node2idx[p] for p in end_prots]
                # only store the length of the shortest path to any of the targets
                any_target_lengths = {}
                for i in range(len(start_idx)):
                    distances = dist_matrix[i][end_idx]
                    min_idx = np.argmin(distances)
                    # convert the distance back to the normalized weights?
                    #min_dist = 10**(-distances[min_idx])
                    min_dist = distances[min_idx]
                    any_target_lengths[(start_prots[i],end_prots[min_idx])] = min_dist
                # now write to file
                write_path_lengths(any_target_lengths, out_file, taxon=taxon, goid=goid, write_type='a')
                 
        else:
            for start_prots, end_prots, out_file in start_end_sets:
                if len(goid_pos) <= 2:
                    out_file = out_file.replace(".txt","%s-%s.txt" % (taxon, goid))
                path_lengths = find_shortest_paths(G, start_prots, end_prots, unweighted=unweighted)
                # only store the length of the shortest path to any of the targets
                any_target_lengths = {}
                old_p1 = ""
                min_p2 = ""
                min_l = 1000
                for (p1,p2), length in sorted(path_lengths.items()):
                    if p1 != old_p1:
                        # store the last iterations min path lengths
                        if min_l != 1000:
                            any_target_lengths[(old_p1,min_p2)] = min_l
                        min_l = 1000
                        old_p1 = p1
                    if float(length) < min_l:
                        min_l = length
                        min_p2 = p2
                # make sure to get the last one
                if min_l != 1000:
                    any_target_lengths[(old_p1,min_p2)] = min_l

                write_path_lengths(any_target_lengths, out_file, taxon=taxon, goid=goid, write_type='a')


def find_shortest_paths_scipy(P, sources, directed=False):
    sparse.csgraph.shortest_path(P, directed=directed, indices=sources)

# first setup the graph



def find_shortest_paths(G, sources, targets, unweighted=True, out_file=None):
    sources = set(sources)
    targets = set(targets)
    print("Finding the shortest paths from %d taxon prots to %d non-taxon prots (%d total)"
          % (len(sources), len(targets), len(sources) * len(targets)))

    path_lengths = {}
    #for i, p1 in enumerate(tqdm(sorted(sources))):
    for i, p1 in enumerate(sorted(sources)):
        for j, p2 in enumerate(sorted(targets)):
            # if the sources and targets are the same, then only consider each pair once
            if len(sources) == len(targets):
                if p2 <= p1:
                    continue
            length = find_shortest_path(G, p1, p2, unweighted=unweighted)
            path_lengths[(p1, p2)] = length
            #print(p1,p2,length)

    if out_file is not None:
        # now write to a file / plot a histogram
        write_path_lengths(path_lengths, out_file)
    return path_lengths


def find_shortest_path(G, source, target, unweighted=True):
    try:
        length = nx.shortest_path_length(
            G, source, target, weight="weight" if unweighted is False else None)
    except nx.exception.NetworkXNoPath:
        length = "inf"
    return length


def write_path_lengths(path_lengths, out_file, taxon='-', goid='-', write_type='w'):
    """
    *write_type*: either 'w' for write, or 'a' for append
    """
    if write_type == 'w':
        print("\tWriting to %s" % (out_file))
    else:
        print("\tAppending to %s" % (out_file))
    utils.checkDir(os.path.dirname(out_file))
    with open(out_file, write_type) as out:
        # lock it to avoid scripts trying to write at the same time
        fcntl.flock(out, fcntl.LOCK_EX)
        out.write(''.join("%s\t%s\t%s\t%s\t%s\n" % (taxon, goid,
            p1, p2, str(length)) for (p1,p2), length in sorted(path_lengths.items())))
        fcntl.flock(out, fcntl.LOCK_UN)


def load_network(net_file=None, version_net_file=None, ss_lambda=None):
    # TODO I should syncronize how the node ids and normalized networks and such are stored 
    if net_file is not None:
        # if the file is already an npz, then just read it
        print("Reading network from %s" % (net_file))
        W = sparse.load_npz(net_file)
        P = alg_utils.normalizeGraphEdgeWeights(W, l=ss_lambda)
        node2idx_file = net_file + "-node-ids.txt"
        print("Reading node names from %s" % (node2idx_file))
        node2idx = {n: int(n2) for n, n2 in utils.readColumns(node2idx_file, 1, 2)}
        idx2node = {n2: n for n, n2 in node2idx.items()}
    else:
        if ss_lambda is None:
            sparse_net_file = version_net_file.replace('.txt', '-normalized.npz')
            node2idx_file = version_net_file.replace(".txt", "-node2idx.txt")
        else:
            sparse_net_file = version_net_file.replace('.txt', '-normalized-l%s.npz' % (
                str(ss_lambda).replace('.', '_')))
            node2idx_file = version_net_file.replace(".txt", "-node2idx-l%s.txt" % (
                str(ss_lambda).replace('.', '_')))
        if os.path.isfile(sparse_net_file) and os.path.isfile(node2idx_file):
            print("Reading network from %s" % (sparse_net_file))
            P = sparse.load_npz(sparse_net_file)
            print("Reading node names from %s" % (node2idx_file))
            node2idx = {n: int(n2) for n, n2 in utils.readColumns(node2idx_file, 1, 2)}
            idx2node = {n2: n for n, n2 in node2idx.items()}
        else:
            print("ERROR: %s network not yet setup. Use run_algs.py.\nQuitting" % (version_net_file))
            sys.exit(1)

    return P, node2idx, idx2node


if __name__ == '__main__':
    # I'm reusing the options from the eval_leave_one_species_out script
    opts = eval_sp.parse_args(sys.argv)
    if opts.taxon is None:
        # TODO use all of the species?
        sys.exit("ERROR: must specify at least one --taxon")

    goterms = alg_utils.select_goterms(
            only_functions_file=opts.only_functions, goterms=opts.goterm) 
    goid_pos, goid_neg = alg_utils.parse_pos_neg_files([opts.pos_neg_file], goterms=goterms) 
    eval_pos, eval_neg = None, None
    if opts.pos_neg_file_eval:
        eval_pos, eval_neg = alg_utils.parse_pos_neg_files([opts.pos_neg_file_eval], goterms=goterms) 

    #version = opts.version[0]
    #alpha = opts.alpha[0]

    for taxon in opts.taxon:
        main(opts.version, opts.exp_name, goid_pos, goid_neg, opts.algorithm,
             taxon=taxon, unweighted=opts.unweighted,
             eval_pos=eval_pos, eval_neg=eval_neg)
