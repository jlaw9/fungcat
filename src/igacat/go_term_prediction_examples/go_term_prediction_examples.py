#!/usr/bin/env 

import os, sys
from optparse import OptionParser
from collections import defaultdict
import networkx as nx
# See the tutorial of obonet here: 
# https://github.com/dhimmel/obonet/blob/master/examples/go-obonet.ipynb
import obonet
import pandas as pd
from tqdm import tqdm


# Guide to GO evidence codes: http://geneontology.org/page/guide-go-evidence-codes
ALL_EVIDENCE_CODES = ['EXP','IDA','IPI','IMP','IGI','IEP','ISS','ISO','ISA','ISM','IGC','IBA','IBD','IKR','IRD','RCA','TAS','NAS','IC','ND','IEA']


def parse_obo_file_and_build_dags(obo_file):
    """
    Parse the GO OBO into a networkx MultiDiGraph using obonet.
    Then construct a DAG for each category using the 'is_a' relationships 
    
    *returns*: a dictionary containing a DAG for each of the 3 GO categories 'C', 'F', and 'P'
    """
    global id_to_name
    global name_to_id
    global goid_to_category  # mapping from a GO term ID and the category it belongs to ('C', 'F' or 'P')

    print("Reading GO OBO file from %s" % (obo_file))
    # obonet returns a networkx MultiDiGraph object containing all of the relationships in the ontology
    graph = obonet.read_obo(obo_file)
    # build a mapping from the GO term IDs to the name of the GO term
    id_to_name = {id_: data['name'] for id_, data in graph.nodes(data=True)}
    name_to_id = {data['name']: id_ for id_, data in graph.nodes(data=True)}
    goid_to_category = {} 
    print("\t%d nodes, %d edges" % (graph.number_of_nodes(),graph.number_of_edges()))

    # make sure this really is a DAG
    if nx.is_directed_acyclic_graph(graph) is False:
        print("\tWarning: graph is not a dag")

    # copied this section from cell 19 of https://github.com/IGACAT/DataPreprocessing/blob/master/scripts/populate_go_terms.ipynb
    # Extract all edges with "is_a" relationship.
    # I did not include "part_of" relationships because the molecular_function and biological_process DAGs are not separate from each other if I do
    is_a_edge_list = []
    for child, parent, key in graph.out_edges(keys=True):
        if key == 'is_a':
            is_a_edge_list.append((child, parent))

    # get a is_a-type edge-induced subgraph 
    is_a_subG = nx.MultiDiGraph(is_a_edge_list)
    full_to_category = {'cellular_component': 'C', 'biological_process': 'P', 'molecular_function' : 'F'}
    go_dags = {}
    # there are 3 weakly_connected_components. One for each category
    for wcc in nx.weakly_connected_components(is_a_subG):
        G = is_a_subG.subgraph(wcc)

        # store this DAG in the dictionary of GO DAGs
        # find the root node 
        root_node = None  # find root_node  (no out_edge) 
        for node in G.nodes():
            if G.out_degree(node) == 0:
                root_node = node
                #print(root_node, id_to_name[node])
                break
        c = full_to_category[id_to_name[root_node]]
        print("\tDAG for %s has %d nodes" % (id_to_name[root_node], len(wcc)))
        go_dags[c] = G

        # also set the category for each GO term
        for n in G.nodes():
            goid_to_category[n] = c

    return go_dags


def parse_gaf_file(gaf_file, pos_neg_ec=[], rem_neg_ec=[], ignore_ec=[]):
    """
    Parse a GAF file containing direct annotations (i.e. annotations have not been propogated up the GO DAG)

    Calls the function setup_evidence_code_categories() to assign each evidence code 
    to either the *pos_neg_ec* set, *rem_neg_ec* set, or the *ignore_ec* set. See that function
    for more description

    Returns:
    *prot_goids_by_c*: for each category ('C', 'F', or 'P'), 
        contains the set of GO term IDs to which each protein is annotated (*pos_neg_ec* codes)
    *goid_prots*: contains the set of proteins annotated to each GO term ID (*pos_neg_ec* codes)
    *goid_rem_neg_prots*: contains the set of proteins annotated to each GO term ID (*rem_neg_ec* codes)
    *all_prots*: all proteins. Used to assign unknowns
    """

    print("Setting up evidence code categories")
    pos_neg_ec, rem_neg_ec, ignore_ec = setup_evidence_code_categories(pos_neg_ec, rem_neg_ec, ignore_ec)

    print("Reading annotations from GAF file %s." % (gaf_file))

    # dictionary with key: uniprotID, val: set of goterm IDs annotated to the protein
    # split by hierarchy/category so we can just pass a given categories's annotations when defining negatives
    prot_goids_by_c = {"C": defaultdict(set), "F": defaultdict(set), "P": defaultdict(set)}
    # dictionary with key: goterm ID, val: set of proteins annotated to the goterm ID 
    goid_prots = defaultdict(set)
    goid_rem_neg_prots = defaultdict(set)
    all_prots = set()
    num_not_ann = 0
    num_pos_neg_ann = 0
    num_rem_neg_ann = 0
    num_ignored_ann = 0 

    # if they pass in a GAF file:
    with open(gaf_file, 'r') as f:
        for line in f:
            cols = line.rstrip().split('\t')
            prot = cols[1]
            all_prots.add(prot)
            goid = cols[4]
            evidence_code = cols[6]
            category = cols[8]
            # for now, ignore cellular component annotations
            if category == "C":
                continue
            # skip NOT annotations for now
            if "NOT" in cols[3]:
                num_not_ann += 1 
                continue

            if evidence_code in ignore_ec:
                num_ignored_ann += 1
            elif evidence_code in pos_neg_ec:
                num_pos_neg_ann += 1
                prot_goids_by_c[category][prot].add(goid)
                goid_prots[goid].add(prot) 
            elif evidence_code in rem_neg_ec:
                num_rem_neg_ann += 1
                goid_rem_neg_prots[goid].add(prot)
            else:
                print("WARNING: evidence_code '%s' not recognized" % (evidence_code))

    print("\t%d NOT annotations ignored" % (num_not_ann)) 
    print("\t%d \"pos_neg_ec\" annotations" % (num_pos_neg_ann))
    print("\t%d \"rem_neg_ec\" annotations" % (num_rem_neg_ann))
    print("\t%d \"ignore_ec\" annotations" % (num_ignored_ann))
    print("\t%d proteins have 1 or more BP annotations" % (len(prot_goids_by_c["P"])))
    print("\t%d proteins have 1 or more MF annotations" % (len(prot_goids_by_c["F"])))

    return prot_goids_by_c, goid_prots, goid_rem_neg_prots, all_prots


def setup_evidence_code_categories(pos_neg_ec=[], rem_neg_ec=[], ignore_ec=[]):
    """
    Assigns each evidence code to either the *pos_neg_ec*, *rem_neg_ec*, or the *ignore_ec* set
    *pos_neg_ec*: a list of GO evidence codes used to assign positive and negative examples.
            If none are specified, all evidence codes not in the two other categories will be put in this category by default.
    *rem_neg_ec*: a list of GO evidence codes used to remove negative examples.
            Specifically, If a protein would be labelled as a negative example for a given term 
            but is annotated with a "rem_neg" evidence code for the term, it is instead labelled as unknown.
            If none are specified, but "pos_neg_ec" codes are given, 
            all codes not in the other two categories will be put in this category by default.
    *ignore_ec*: a list of GO evidence codes to ignore completely when parsing the GAF file.
            If both --pos-neg-ec and --rem-neg-ec codes are given, everything else will be ignored by default.
            ND is always ignored.

    *returns*: *pos_neg_ec*, *rem_neg_ec*, *ignore_ec* 
    """
    # the ND annotation means there is no data available for this protein. 
    # more information about the ND annotation is available here: http://geneontology.org/page/nd-no-biological-data-available
    if "ND" not in ignore_ec:
        print("\tIngoring the evidence code 'ND' because it means there is no data available for this protein")
        ignore_ec.append("ND")

    # set the positive codes to all of them by default
    # use lists instead of sets here to keep the original order
    if len(pos_neg_ec) == 0:
        # don't use sets to keep the order of the codes
        pos_neg_ec = [c for c in ALL_EVIDENCE_CODES
                      if c not in rem_neg_ec and
                      c not in ignore_ec]
        #pos_neg_ec = set(ALL_EVIDENCE_CODES).difference(set(rem_neg_ec)) \
        #                                                 .difference(set(ignore_ec))
    # if 1 or more positive evidence codes are given, but no non-negative codes are given,
    # set the rest of the codes to be non-negative by default
    elif len(rem_neg_ec) == 0:
        rem_neg_ec = [c for c in ALL_EVIDENCE_CODES
                      if c not in pos_neg_ec and
                      c not in ignore_ec]
    # if 1 or more positive and 1 or more non-negative codes are given,
    # set the rest to be ignored by default
    else:
        ignore_ec = [c for c in ALL_EVIDENCE_CODES
                     if c not in pos_neg_ec and
                     c not in rem_neg_ec]

    print()
    print("pos_neg_ec (used to assign positive and negative examples):" +
          "\n\t'%s'" % ("','".join(pos_neg_ec))) 
    print("rem_neg_ec (used to remove negative examples):" +
          "\n\t'%s'" % ("','".join(rem_neg_ec))) 
    print("ignore_ec (ignored completely when assigning examples):" +
          "\n\t'%s'" % ("','".join(ignore_ec)))
    print()

    # make sure the sets are non-overlapping
    if len(set(pos_neg_ec).intersection(set(rem_neg_ec))) != 0 or \
       len(set(pos_neg_ec).intersection(set(ignore_ec))) != 0 or \
       len(set(rem_neg_ec).intersection(set(ignore_ec))) != 0:
        sys.stderr.write("ERROR: the three sets are not disjoint. " +
                         "Please ensure the three input sets have no overlapping evidence codes.\n")
        sys.exit(1)

    return pos_neg_ec, rem_neg_ec, ignore_ec


def extract_high_freq_goterms(G, goids, annotated_prots, cutoff=1000):
    """
    *G*: GO DAG (networkx DiGraph) with prot->goid edges for each protein's annotations
    returns a set of GO terms with > cutoff proteins annotated to it 
    """
    high_freq_go_terms = set() 
    for goid in tqdm(goids):
        anc = nx.ancestors(G, goid)
        # the number of positive annotations for this GO term is the number of proteins that can reach this GO term ID in the gene-goid graph
        # meaning the number of proteins annotated to this term plus those annotated to an ancestral, more specific term
        if len(anc.intersection(annotated_prots)) > cutoff:
            high_freq_go_terms.add(goid) 

    return high_freq_go_terms


def build_gene_goterm_graph(go_dag, goid_prots):
    """
    For every protein, add an edge from the protein to the GO term IDs to which it's annotated
    *go_dag*: networkx DiGraph DAG containing the is_a edges in the GO DAG 
    *goid_prots*: contains the set of proteins annotated to each GO term ID

    *returns*: the resulting gene-goterm graph (networkx DiGraph), and the graph reversed.
    """

    G = nx.DiGraph()
    G.add_edges_from(go_dag.edges())

    # revG is a copy of the annotation graph G with the GO DAG reversed
    revG = nx.reverse(G, copy=True)

    # set all of the current nodes as goids
    #nx.set_node_attributes(G, 'goid', 'type')

    # For every GO term ID, add an edge in the graph from the proteins annotated to the GO term ID, to the GO term ID
    # This graph allows us to get all of the proteins annotated to descendants (more specific terms) of a term
    for goid in go_dag.nodes():
        for prot in goid_prots[goid]:
            # add an edge from the protein to the GO term its annotated to
            G.add_edge(prot, goid)
            revG.add_edge(prot, goid)

    print("\t%d nodes, %d edges" % (G.number_of_nodes(),G.number_of_edges()))

    return G, revG


def assign_pos_neg(goid, G, revG, annotated_prots, rem_negG=None):
    """
    This function assigns the set of positive and negative proteins for a given GO term ID.
    Specifically, for the given GO term t, we define a gene/protein g as a 
    - positive if g is directly annotated to t or to a descendant of t (more specific term) in the GO DAG
    - negative if g is not annotated to t or an ancestor or descendant of t in the GO DAG, but also has at least 1 other annotation
    - unknown if g is neither a positive nor a negative meaning it has no annotations, 
      or is annotated to an ancestor of t (more general term) in the GO DAG

    Parameters:
    *goid*: GO term for which to assign positives and negatives
    *G*: GO DAG with prot->goid edges for each protein's annotations
    *revG*: reverse of G. Used to find all of the proteins annotated to descendant or less-specific GO terms 
    *annotated_prots*: all proteins with at least one direct annotation (in the GO category of the given GO term). 
        Used to assign negatives and get the protein nodes from G and revG
    *rem_negG*: version of the annotation graph G which contains the rem_neg_ec annotations to remove negative examples

    Returns:
    *positives*: the set of proteins labelled as positives
    *negatives*: the set of proteins labelled as negatives 
    """

    # positives are all of the proteins can reach this GO term.
    positives = set(nx.ancestors(G, goid)).intersection(annotated_prots)
    # proteins that can be reached from this term are unknowns 
    unknowns = set(nx.ancestors(revG, goid)).intersection(annotated_prots)
    ## if this node is directly annotated to the term, it's a positive
    #unknowns.difference_update(positives)

    if rem_negG is not None:
        # if the protein is annotated to the term, or a more specific term, with a non-negative (remove negative) evidence code,
        # don't use it as a negative
        rem_negs = set(nx.ancestors(rem_negG, goid)).intersection(annotated_prots)
        unknowns.update(rem_negs)

    # negatives are all of the proteins with an annotation that are not an ancestor, or descendant 
    negatives = annotated_prots.difference(positives) \
                               .difference(unknowns)

    return positives, negatives


def assign_all_pos_neg(high_freq_goids, G, revG, annotated_prots, all_prots, rem_negG=None, verbose=False):
    """
    Assigns each gene as a positive/negative/unknown example for each GO term. 

    Parameters: 
    *high_freq_goids*: goids for which to get positives and negatives. Should all belong to a single category
    *G*: annotation graph
    *revG*: annotation graph with GO DAG reversed
    *annotated_prots*: all proteins with a direct annotation (in the GO category of the high_freq_goids). Used to assign negatives
    *all_prots*: all proteins. Used to assign unknowns
    *rem_negG*: version of the annotation graph G which contains the rem_neg_ec annotations to remove negative examples
    *verbose*: print the # of positives, negatives and unknowns for each GO term

    Returns:
    *goid_pos*: dictionary of a set of positive examples for each GO term
    *goid_neg*: dictionary of a set of negative examples for each GO term
    *goid_unk*: dictionary of a set of unknown examples for each GO term
    """
    global id_to_name, name_to_id

    print("Getting positives and negatives for %d GO terms" % (len(high_freq_goids)))

    # dictionaries containing the set of positives and negatives respectively for each GO term ID
    goid_pos = {}
    goid_neg = {}
    goid_unk = {}
    # for each GO term, get the set of positives and the set of negatives, and store them in a dictionary
    for goid in tqdm(sorted(high_freq_goids)):
        positives, negatives = assign_pos_neg(goid, G, revG, annotated_prots, rem_negG=rem_negG)
        goid_pos[goid] = positives
        goid_neg[goid] = negatives
        goid_unk[goid] = all_prots.difference(positives).difference(negatives) 
        if verbose is True:
            tqdm.write("\t%d positives, %d negatives, %d unknowns for %s (%s)" % (len(positives), len(negatives), len(goid_unk[goid]), id_to_name[goid], goid))

    return goid_pos, goid_neg, goid_unk


def build_pos_neg_table(high_freq_goids, goid_pos, goid_neg, goid_unk, summary_only=False):
    """
    Builds a table with a positive/negative/unknown (1/-1/0) assignment for each gene-GO term pair. 
    Rows are the genes and columns are the given high_freq_goids (GO terms with > cutoff proteins annotated) 

    Parameters: 
    *high_freq_goids*: goids for which to get positives and negatives. Should all belong to a single category
    *goid_pos*: positive examples for each GO term
    *goid_neg*: negative examples for each GO term
    *goid_unk*: unknown examples for each GO term
    *summary_only*: build and return only the summary table

    Returns:
    *df*: the table as a pandas DataFrame 
    *df_summary*: a table containing the # of positive, negative and unknown examples for each GO term
    """
    global id_to_name, name_to_id, goid_to_category

    if summary_only is False:
        print("Building a table with positive/negative/unknown assignments for each protein-goterm pair")
        # build a table with the first column being the genes, and a column for each of the terms with > cutoff annotations indicating 1/-1/0 assignment for each gene
        pos_neg_table = defaultdict(dict)
        # build a double dictionary with either 1, -1 or 0 for each GO term protein pair
        # TODO there must be a better pandas method to construct the table
        for goid in tqdm(high_freq_goids):
            for prot in goid_pos[goid]:
                pos_neg_table[goid][prot] = 1
            for prot in goid_neg[goid]:
                pos_neg_table[goid][prot] = -1
            # unknowns are everything that is not a positive or negative
            for prot in goid_unk[goid]:
                pos_neg_table[goid][prot] = 0

        df = pd.DataFrame(pos_neg_table)

    df_summary = pd.DataFrame({
        "GO term name": {goid: id_to_name[goid] for goid in high_freq_goids},
        "GO category": {goid: goid_to_category[goid] for goid in high_freq_goids},
        "# positive examples": {goid: len(pos) for goid, pos in goid_pos.items()}, 
        "# negative examples": {goid: len(neg) for goid, neg in goid_neg.items()}, 
        "# unknown examples": {goid: len(unk) for goid, unk in goid_unk.items()}
        })
    # set the order of the columns
    cols = ["GO term name", "GO category", "# positive examples", "# negative examples", "# unknown examples"]
    df_summary = df_summary[cols] 
    df_summary.index.rename("GO term", inplace=True)

    if summary_only is False:
        return df, df_summary
    else:
        return df_summary


def main(obo_file, gaf_file, out_pref, cutoff=1000, write_table=False,
        pos_neg_ec=[], rem_neg_ec=[], ignore_ec=[]):
    # first parse the gaf and obo files
    direct_prot_goids_by_c, direct_goid_prots, direct_goid_rem_neg_prots, all_prots = parse_gaf_file(
            gaf_file, pos_neg_ec, rem_neg_ec, ignore_ec)
    go_dags = parse_obo_file_and_build_dags(obo_file)

    # keep track of the summary stats for each category, and combine them into one table in the end
    df_summaries = pd.DataFrame()

    # assign the positives, negatives and unknowns for biological process and molecular function
    for c in ["P", "F"]:
        print("Category: %s" % (c))
        print("Building the gene-goterm graph")
        G, revG = build_gene_goterm_graph(go_dags[c], direct_goid_prots)
        rem_negG = None
        if len(direct_goid_rem_neg_prots) > 0:
            # the remove-negative annotations also need to be be propagated, so build an annotation graph for them here
            rem_negG, rem_neg_revG = build_gene_goterm_graph(go_dags[c], direct_goid_rem_neg_prots)
        #print("# of prots with at least 1 %s annotation: %d" % (c, len(prot_goids)))
        #print("# of %s GO terms with at 1 protein annotated to it: %d" % (c, len(goid_prots)))

        print("Extracting GO terms with > %d annotations" % (cutoff))
        annotated_prots = set(direct_prot_goids_by_c[c].keys())
        high_freq_goids = extract_high_freq_goterms(G, go_dags[c].nodes(), annotated_prots, cutoff=cutoff)
        # also remove biological process, cellular component and molecular function
        high_freq_goids.difference_update(set([name_to_id[name] for name in ["cellular_component", "biological_process", "molecular_function"]]))
        print("\t%d (out of %d) GO terms have > %d proteins annotated to them" % (len(high_freq_goids), go_dags[c].number_of_nodes(), cutoff))

        # keep track of the set of proteins with at least 1 annotation in this category to assign negatives later
        goid_pos, goid_neg, goid_unk = assign_all_pos_neg(high_freq_goids, G, revG, annotated_prots, all_prots, rem_negG=rem_negG)

        # now write it to a file
        category = {"C": "cc", "P": "bp", "F": "mf"}
        if write_table is True:
            # build a table containing a positive/negative/unknown assignment for each protein-goterm pair
            df, df_summary = build_pos_neg_table(high_freq_goids, goid_pos, goid_neg, goid_unk)
            # combine the summary stats for all categories into one table
            df_summaries = pd.concat([df_summaries, df_summary])

            out_file = "%spos-neg-%s-%d.tsv" % (out_pref, category[c], cutoff)
            print("Writing table containing positive/negative/unknown assignments to %s" % (out_file))
            df.to_csv(out_file, sep="\t")
        else:
            # build a summary table of the pos/neg/unk assignments
            df_summary = build_pos_neg_table(high_freq_goids, goid_pos, goid_neg, goid_unk, summary_only=True)
            # combine the summary stats for all categories into one table
            df_summaries = pd.concat([df_summaries, df_summary])
            out_file = "%spos-neg-%s-%d-list.tsv" % (out_pref, category[c], cutoff)
            print("Writing file containing positive/negative assignments to %s" % (out_file))
            with open(out_file, 'w') as out:
                out.write("#goid\tpos/neg assignment\tprots\n")
                for goid in high_freq_goids:
                    out.write("%s\t1\t%s\n" % (goid, ','.join(goid_pos[goid])))
                    out.write("%s\t-1\t%s\n" % (goid, ','.join(goid_neg[goid])))

    output_summary_file = "%spos-neg-%d-summary-stats.tsv" % (out_pref, cutoff)
    # maybe make this into an option later instead of always writing it
    #if output_summary_file is not None:
    print("Writing summary table of # of positive, negative and unknown examples for each GO term to: %s" % (output_summary_file))
    df_summaries.to_csv(output_summary_file, sep='\t')


def parse_args(args):
    ## Parse command line args.
    description = """
This script takes the annotations in a GAF file, and the GO DAG and assigns 
every gene as either a positive (1), negative (-1) or unknown (0) for each GO term with > cutoff annotations.
Writes two tab-separated tables containing the assignments, one for BP and one for MF, where the rows are genes, 
and the columns are GO term IDs. Also writes a summary statistics table
"""
    usage = '%prog [options] '
    parser = OptionParser(usage=usage, description=description)
    parser.add_option('-g', '--gaf-file', type='string',
                      help="File containing GO annotations in GAF format. Required")
    parser.add_option('-b', '--obo-file', type='string', 
                      help="GO OBO file which contains the GO DAG. Required")
    #parser.add_option('-n', '--negatives', type='string', default='non-ancestral',
    #                  help="Types of negatives to generate. Options are: '%s'. Default = 'non-ancestral', See the README file for descriptions of these options." % ("', '".join(NEGATIVES_OPTIONS)))
    parser.add_option('-c', '--cutoff', type='int', default=1000,
                      help="GO terms having > cutoff positive instances (proteins) are kept. Default=1000")
    parser.add_option('-o', '--out-pref', type='string', 
                      help="Prefix used to write a table of positives, negatives, and unknowns for each GO category." +
                      "Writes an output file for BP and MF: <out-pref>pos-neg-<cutoff>-P.tsv and <out-pref>pos-neg-<cutoff>-F.tsv")
    # writing the big pos/neg/unk assignment matrix is taking too long. 
    # instead, write the pos/neg prots for each GO term to a file
    parser.add_option('', '--write-table', action='store_true', default=False,
                      help="write the pos/neg/unk assignments to a table rather than the default comma-separated list of prots")
    parser.add_option('', '--pos-neg-ec', type='string',
                      help="Comma-separated list of evidence codes used to assign positive and negative examples. " +
                      "If none are specified, all codes not in the two other categories " + 
                      "(--rem-neg-ec and --ignore-ec) will be used by default.")
    parser.add_option('', '--rem-neg-ec', type='string',
                      help="Comma-separated list of evidence codes used to remove negative examples. " + 
                      "Specifically, If a protein would be labelled as a negative example for a given term " + 
                      "but is annotated with a 'rem_neg' evidence code for the term, it is instead labelled as unknown. " +
                      "If none are specified, but --pos-neg-ec codes are given, " +
                      "all codes not in the other two categories will be put in this category by default.")
    parser.add_option('', '--ignore-ec', type='string',
                      help="Comma-separated list of evidence codes where annotations with the specified codes will be ignored when parsing the GAF file. " +
                      "For example, specifying 'IEA' will skip all annotations with an evidence code 'IEA'. " +
                      "If both --pos-neg-ec and --rem-neg-ec codes are given, everything else will be ignored by default.")

    (opts, args) = parser.parse_args(args)

    if opts.gaf_file is None or opts.obo_file is None or opts.out_pref is None:
        parser.print_help()
        sys.exit("\n--gaf-file (-g), --obo-file (-b), and --out-pref (-o) are required")

    # make sure all of the specified codes are actually GO evidence codes
    codes = []
    for codes_option in [opts.pos_neg_ec, opts.rem_neg_ec, opts.ignore_ec]:
        if codes_option is not None:
            codes += codes_option.split(',')
    non_evidence_codes = set(codes).difference(set(ALL_EVIDENCE_CODES))
    if len(non_evidence_codes) > 0:
        sys.stderr.write("ERROR: the specified code(s) are not GO evidence codes: '%s'\n" % ("', '".join(non_evidence_codes)))
        sys.stderr.write("Accepted evidence codes: '%s'\n" % ("', '".join(ALL_EVIDENCE_CODES)))
        sys.exit(1)

    # check if the output prefix is writeable
    out_dir = os.path.dirname(opts.out_pref)
    if not os.path.isdir(out_dir):
        sys.stderr.write("ERROR: output directory %s specified by --out-pref doesn't exist\n" % (out_dir))
        sys.exit(1)

    return opts, args


if __name__ == "__main__":
    print("Running %s" % (' '.join(sys.argv)))
    opts, args = parse_args(sys.argv)
    pos_neg_ec = [] if opts.pos_neg_ec is None else opts.pos_neg_ec.split(',') 
    rem_neg_ec = [] if opts.rem_neg_ec is None else opts.rem_neg_ec.split(',') 
    ignore_ec = [] if opts.ignore_ec is None else opts.ignore_ec.split(',') 
    main(opts.obo_file, opts.gaf_file, opts.out_pref, cutoff=opts.cutoff, write_table=opts.write_table,
         pos_neg_ec=pos_neg_ec, rem_neg_ec=rem_neg_ec, ignore_ec=ignore_ec)
