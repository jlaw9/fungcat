#!/usr/bin/python

# code to post toxcast results to graphspace

#print("Importing Libraries"))

from collections import defaultdict
from optparse import OptionParser
from graphspace_python.api.client import GraphSpace
from graphspace_python.graphs.classes.gsgraph import GSGraph
# GSGraph already implements networkx
import networkx as nx
import utils.file_utils as utils
import fungcat_settings as f_settings
# for string networks
import src.algorithms.setup_sparse_networks as setup
import os, sys
import utils.graphspace.post_to_graphspace as gs
import algorithms.gain_scipy.alg_utils as alg_utils
from pandas import read_csv


evidence_code_name = {
    "EXP": "Inferred from Experiment",
    "IDA": "Inferred from Direct Assay",
    "IPI": "Inferred from Physical Interaction",
    "IMP": "Inferred from Mutant Phenotype",
    "IGI": "Inferred from Genetic Interaction",
    "IEP": "Inferred from Expression Pattern",
    "ISS": "Inferred from Sequence or structural Similarity",
    "ISO": "Inferred from Sequence Orthology",
    "ISA": "Inferred from Sequence Alignment",
    "ISM": "Inferred from Sequence Model",
    "IGC": "Inferred from Genomic Context",
    "IBA": "Inferred from Biological aspect of Ancestor",
    "IBD": "Inferred from Biological aspect of Descendant",
    "IKR": "Inferred from Key Residues",
    "IRD": "Inferred from Rapid Divergence",
    "RCA": "Inferred from Reviewed Computational Analysis",
    "TAS": "Traceable Author Statement",
    "NAS": "Non-traceable Author Statement",
    "IC" : "Inferred by Curator",
    "ND" : "No biological Data available",
    "IEA": "Inferred from Electronic Annotation",
}

evidence_code_type = {"IEA": "electronic"}
for code in ["EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"]:
    evidence_code_type[code] = "experimental"
# Author Statement evidence codes are grouped in here as well
for code in ["ISS", "ISO", "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA", "NAS"]:
    evidence_code_type[code] = "computational"

# color nodes according to their type
node_type_color = {
    "prediction": "#d88c00",  # orange
    "annotation": "#40aff9",  # blue
    # also keep track of the negatives
    "neg-annotation": "#8f68a5",  # dark purple
    "non-taxon-annotation": "#8ec67b",  # green
    "non-taxon-neg-annotation": "#54575b",  # default grey 
    "default": "#D8D8D8",  # default grey - background-color
}
annotation_type_styles = {
    # all experimental evidence codes. Maroon double border, square
    "experimental": {"border-style": "double", "border-width": "10", "border-color": "#cc5800", "shape": "rectangle"},  
    # computational analysis evidence codes (still with some curator input). Black double border, square
    "computational": {"border-style": "double", "border-width": "10", "border-color": "#3a3835", "shape": "rectangle"},  
    # Inferred by Electronic Annotation. For now, just use the defaults
    "electronic": {"shape": "rectangle"},  
    #"electronic": {},  
}
edge_type_color = {
    "default": "#6b6b6b",
    "string": "#24bf1c",
}

NAME_TO_SHORTNAME = {
    "Neisseria gonorrhoeae FA 1090"                             : "Ng",
    "Peptoclostridium difficile / Clostridioides difficile 630" : "Cd",
    "Helicobacter pylori 85962"                                 : "Hp",
    "Klebsiella pneumoniae"                                     : "Kp",
    "Enterococcus faecium DO"                                   : "Ef",
    "Escherichia coli K-12"                                     : "Ec",
    "Haemophilus influenzae RD KW20"                            : "Hi",
    "Bordetella pertussis Tohama I"                             : "Bp",
    "Burkholderia cepacia"                                      : "Bc",
    "Clostridium botulinum A str. Hall"                         : "Cb",
    "Acinetobacter baumannii"                                   : "Ab",
    "Staphylococcus aureus"                                     : "Sa",
    "Vibrio cholerae O1 biovar El Tor str. N16961"              : "Bc",
    "Yersinia pestis"                                           : "Yp",
    "Streptococcus pyogenes"                                    : "Sp",
    "Pseudomonas aeruginosa"                                    : "Pa",
    "Salmonella typhimurium / Salmonella enterica"              : "Se",
    "Shigella dysenteriae serotype 1 (strain Sd197)"            : "Sd",
    "Mycobacterium tuberculosis"                                : "Mt",
}

alg_names = {
    "sinksource": "SinkSource",
    "localplus": "Local+",
    "sinksourceplus": "SinkSource+",
    "genemania": "GeneMANIA",
    "birgrank": "BirgRank",
}


def get_mat_neighbors(A, n):
    neighbors = A[n].nonzero()[1]
    return neighbors

def setup_post_to_graphspace(
        version, algorithm, exp_name, selected_goid, postfix='', tags=None,
        pos_neg_files=None, taxon=None, goid_summary_file=None,
        num_neighbors=1, nodes_to_post=None,
        pos_neg_files_eval=None, keep_ann=False,
        alpha=1, eps=0.0001, maxi=20):
    INPUTSPREFIX, RESULTSPREFIX, NETWORK, selected_strains = f_settings.set_version(version)

    #selected_goid = "15643"  # toxic substance binding
    #selected_goid = "9405"  # pathogenesis
    #selected_goid = "98754"  # detoxification
    selected_goname = None
    # build a dictionary of the evidencecode for each prot
    uniprot_to_evidencecode = defaultdict(set)
    annotated_prots = set()
    neg_prots = set()
    if goid_summary_file is not None:
        df_summary = read_csv(goid_summary_file, sep='\t')
        goid_names = dict(zip(df_summary['GO term'], df_summary['GO term name']))
        #goid_num_anno = dict(zip(df_summary['GO term'], df_summary['# positive examples']))
        print("GO name: %s" % (goid_names[selected_goid]))
        selected_goname = goid_names[selected_goid].replace(' ','-')[0:20]
    #goid_to_goname = {}
    if pos_neg_files is not None:
        goid_pos, goid_neg = alg_utils.parse_pos_neg_files(pos_neg_files, goterms=set([selected_goid])) 
        annotated_prots = goid_pos[selected_goid]
        neg_prots = goid_neg[selected_goid]
    if pos_neg_files_eval is not None:
        goid_pos, goid_neg = alg_utils.parse_pos_neg_files(pos_neg_files_eval, goterms=set([selected_goid])) 
        eval_annotated_prots = goid_pos[selected_goid]
        eval_neg_prots = goid_neg[selected_goid]
    # load the GAIN propagation to get the evidence code 
    for orf, goid, goname, hierarchy, evidencecode, annotation_type in utils.readColumns(f_settings.GOA_ALL_FUN_FILE_NOPARTOF, 1,2,3,4,5,6):
        if selected_goid[:3] == "GO:":
            goid = "GO:" + "0"*(7-len(goid)) + goid
        if goid != selected_goid:
            continue
        selected_goname = goname.replace(' ','-')[0:20]
        # TODO figure out how to represent negatives 
        if annotation_type != '1':
            continue 

        if pos_neg_files is None:
            annotated_prots.add(orf) 
        uniprot_to_evidencecode[orf].add(evidencecode)
        #goid_to_goname[goid] = goname
    # limit it to the current taxon
    if taxon is not None:
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
        # also limit the proteins to those in the network
        print("\t%d prots for taxon %s." % (len(taxon_prots), taxon))
        if keep_ann:
            non_taxon_annotated_prots = set(annotated_prots)
            non_taxon_neg_prots = set(neg_prots)
        else:
            non_taxon_annotated_prots = set(annotated_prots) - taxon_prots 
            non_taxon_neg_prots = set(neg_prots) - taxon_prots
        print("\t%d non-taxon pos, %d non-taxon neg" % (len(non_taxon_annotated_prots), len(non_taxon_neg_prots)))
        if pos_neg_files_eval is not None:
            print("Using eval for taxon positives and negatives")
            annotated_prots = set(eval_annotated_prots) & taxon_prots 
            neg_prots = set(eval_neg_prots) & taxon_prots
        else:
            annotated_prots = set(annotated_prots) & taxon_prots 
            neg_prots = set(neg_prots) & taxon_prots
        print("\t%d taxon pos, %d taxon neg" % (len(annotated_prots), len(neg_prots)))

    print("\t%d annotated prots for %s (%s)" % (len(annotated_prots), selected_goname, selected_goid))


    #print("Reading network from %s" % (NETWORK))
    #G = nx.Graph(name=version)
    #G = nx.read_weighted_edgelist(NETWORK, delimiter="\t", create_using=G)
    #print("\t%d nodes, %d edges" % (G.number_of_nodes(), G.number_of_edges()))
    out_pref_net = "%s/sparse-nets/" % (INPUTSPREFIX)
    utils.checkDir(out_pref_net)
    # build the file containing the sparse networks
    # TODO make options for this
    sparse_nets, net_names, prots = setup.create_sparse_net_file(
        version, out_pref_net, selected_strains=selected_strains,
        string_nets=setup.CORE_STRING_NETWORKS if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] else [], 
        string_cutoff=f_settings.VERSION_STRING_CUTOFF[opts.version] if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version] else 400,
        forcenet=False)

    # load the network
    if opts.net_file is not None:
    #    net_weights = [4.082354583947708, 0.6336420366998545, 0.2946449160659123, 2.909570093171721, 0.6136365533078896, 0.2834763821632677]
        print("\tloading network from %s" % (opts.net_file))
        W, prots = alg_utils.setup_sparse_network(opts.net_file)
        # also check for network weights
        net_weight_file = opts.net_file + "-net-weights.txt"
        if os.path.isfile(net_weight_file):
            print("\tloading network weights from %s" % (net_weight_file)) 
            weights = {net: float(w) for net, w in utils.readColumns(net_weight_file, 1, 2)}
            print("\t%s" % (', '.join(["%s: %0.3f" % (net, weights[net]) for net in net_names if net in weights])))
            net_weights = []
            indices = []
            for i, net in enumerate(net_names):
                if net in weights:
                    indices.append(i) 
                    net_weights.append(weights[net])
            #net_weights = [float(w) for w in utils.readItemList(net_weight_file, 2)]
    #if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version]:
    #    print("combining networks using hard-coded weights")
    #    indices = [0, 1, 2, 3, 5, 6]
    #    # network weights: 2018_06-seq-sim-e0_1-net.txt: 4.082354583947708, neighborhood: 0.6336420366998545, fusion: 0.2946449160659123, cooccurence: 2.909570093171721, experiments: 0.6136365533078896, database: 0.2834763821632677
    #    # TODO these weights are hard-coded for the GO term SOS response
    #    net_weights = [4.082354583947708, 0.6336420366998545, 0.2946449160659123, 2.909570093171721, 0.6136365533078896, 0.2834763821632677]
    #    combined_network = net_weights[0]*sparse_nets[indices[0]]
    #    for i in range(1,len(indices)):
    #        combined_network += net_weights[i]*sparse_nets[indices[i]] 
    #    W = combined_network
    else:
        W = sparse_nets[0]
    #print("converting to networkx graph (TODO just use sparse matrix)")
    #G = nx.convert_matrix.from_scipy_sparse_matrix(combined_network)
    #G = nx.from_scipy_sparse_matrix(combind_network)
    net_nodes = set(prots)

    #undirG = G.to_undirected()
    annotated_prots = annotated_prots & set(net_nodes)
    neg_prots = neg_prots & set(net_nodes)

    print("\t%d annotated prots are in the network" % (len(annotated_prots)))

    if algorithm not in f_settings.ALGORITHM_OPTIONS:
        predictions_file = "%s/all/%s/%s/pred-a%s-eps%s-maxi%s.txt" % (
            RESULTSPREFIX, algorithm, exp_name,
            str(alpha).replace('.','_'), str(eps).replace('.','_'), str(maxi))
        #if not os.path.isfile(predictions_file):
        #    predictions_file = "%s/all/%s/%s/pred-a1_0-eps0_0001.txt" % (RESULTSPREFIX, algorithm, exp_name)
        #    print("checking %s" % (predictions_file))
        if not os.path.isfile(predictions_file):
            print("%s does not exist." % (predictions_file))
            predictions_file = "%s/all/%s/%s/pred-results.txt" % (RESULTSPREFIX, algorithm, exp_name)
            print("\tchecking %s" % (predictions_file))
    else:
        predictions_file = "%s/all/%s/%s/db-%s-predictions.txt" % (RESULTSPREFIX, algorithm, exp_name, exp_name)
    # TODO GAIN provides a local confience score and a cut confidence score. How should I decide what to use as a cutoff for predictions?
    #conf_cutoff = 0.2
    conf_cutoff = -1
    predicted_prots = set() 
    pred_local_conf = {}
    # TODO fix the cut conf
    pred_cut_conf = {}
    ranks = {}
    curr_rank = 0
    first_zero_rank = None
    print("Reading predictions from %s" % (predictions_file))
    for rank, (gene, goid, local_conf) in enumerate(utils.readColumns(predictions_file, 1,2,3)):
        if algorithm not in f_settings.ALGORITHM_OPTIONS:
            # swap the gene and goid column
            goid, gene = gene, goid
        if goid != selected_goid:
            continue
        # TODO this shouldn't happen
        #if not G.has_node(gene):
        if gene not in net_nodes:
            print("%s not in network" % (gene))
            sys.exit()
            continue
        if taxon is not None:
            if gene in taxon_prots:
                curr_rank += 1 
                ranks[gene] = curr_rank
                if float(local_conf) == 0 and first_zero_rank is None:
                    first_zero_rank = curr_rank 
        else:
            ranks[gene] = rank

        if float(local_conf) >= conf_cutoff \
                and gene not in non_taxon_annotated_prots \
                and gene not in non_taxon_neg_prots:
            predicted_prots.add(gene)
            # move the score between 0 and 1 if it's genemania (normally between -1 and 1)
            # as the score is used to set the opacity
            if algorithm == "genemania":
                pred_cut_conf[gene] = local_conf 
                local_conf = ((float(local_conf) - -1) / float(1--1)) * (1-0) + 0
            pred_local_conf[gene] = local_conf 

    if len(predicted_prots) == 0:
        print("Warning: Scoress for GO term %s not found in scores/predictions file %s" % (selected_goid, predictions_file))
    else:
        print("\t%d predicted prots" % (len(predicted_prots)))
    print("Rank of first zero score: %d" % (first_zero_rank))
    print("Ranks of left-out positives:")
    for gene in sorted(annotated_prots, key=ranks.get):
        print("%s\t%d" % (gene, ranks[gene]))
    print("Including top 30 ranked-proteins of left-out species")
    top_30 = set(sorted(taxon_prots & set(ranks.keys()), key=ranks.get)[:30])
    #for gene in sorted(top_30):
    #    print("%s\t%d" % (gene, ranks[gene]))

    if taxon is not None:
        print("Getting the induced subgraph of the neighbors of the %d annotated nodes" % (len(annotated_prots)))
        prededges = set()
        if nodes_to_post is not None: 
            print("Getting neighbors of %s" % (', '.join(nodes_to_post)))
            nodes_to_add_neighbors = set(nodes_to_post) 
        else:
            nodes_to_add_neighbors = annotated_prots.copy() | top_30
        node2idx = {n: i for i, n in enumerate(prots)}
        for i in range(opts.num_neighbors):
            #print("Adding neighbors %d" % (i+1))
            curr_nodes_to_add_neighbors = nodes_to_add_neighbors.copy()
            nodes_to_add_neighbors = set() 
            print("adding %sneighbors of %d nodes" % ("positive ", len(curr_nodes_to_add_neighbors)))
            for u in curr_nodes_to_add_neighbors:
                #neighbors = set(nx.all_neighbors(G, u))
                neighbors = set([prots[v] for v in get_mat_neighbors(W, node2idx[u])])
                if opts.node_to_post is None:
                    # UPDATE 2018-10: try adding just the positive neighbors of the node
                    # TODO make this a command-line option
                    neighbors = neighbors & (non_taxon_annotated_prots | annotated_prots | top_30)
                #if len(neighbors) > 15 and nodes_to_post is None:
                #    print("\tskipping adding neighbors of %s. len(neighbors): %d" % (u, len(neighbors)))
                #    continue
                nodes_to_add_neighbors.update(neighbors)
                prededges.update(set([(u,v) for v in neighbors]))
    else:
        print("Getting the induced subgraph of the %d annotated and %d predicted proteins" % (len(annotated_prots), len(predicted_prots)))
        print("not yet implemented. quitting")
        sys.exit()
    #    prededges = set(G.subgraph(annotated_prots.union(predicted_prots)).edges())
    prededges = set([tuple(sorted((u,v))) for u,v in prededges])
    # TODO I should also show the disconnected nodes
    prednodes = set([n for edge in prededges for n in edge])

    print("\t%d nodes, %d edges" % (len(prednodes), len(prededges)))
    if len(prededges) > 1000 or len(prednodes) > 500:
        print("\nToo many nodes/edges. Not posting to GraphSpace. Quitting")
        sys.exit()

    #graph_attr_file = ""
    #graph_attr, attr_desc = readGraphAttr()
    # add the edge weight from the network to attr_desc which will be used for the popup
    # set the edges as the neighbors of the annotated genes
    #prededges = set()
    # get the induced subgraph of the annotated nodes and predicted nodes
    #for n in func_prots:
    #    if not G.has_node(n):
    #        continue
    #    for neighbor in G.neighbors(n):
    #        prededges.add((n, neighbor))

    graph_attr = {n: {} for n in prednodes}
    attr_desc = {n: {} for n in prednodes}

    prot_species_file = f_settings.VERSION_UNIPROT_TO_SPECIES[version]
    print("Reading gene names and species for each protein from %s" % (prot_species_file))
    prot_species = utils.readDict(prot_species_file, 1, 2)
    uniprot_to_gene = utils.readDict(prot_species_file, 1, 4)
    # there can be multiple gene names. Just show the first one for now
    uniprot_to_gene = {n:gene.split(' ')[0] for n,gene in uniprot_to_gene.items()}
    node_labels = {} 

    # for each node, add the prediction values
    for n in prednodes:
        # set the name of the node to be the gene name and add the k to the label
        gene_name = uniprot_to_gene.get(n,n)
        species_short_name = NAME_TO_SHORTNAME[f_settings.TAX_TO_NAME[prot_species[n]]]
        # add the species to the front of the gene name
        label = "%s-%s" % (species_short_name, gene_name)
        uniprot_to_gene[n] = label
        #node_labels[n] = "%s\n%d" % (label, min(ranks[n], 43)) if n in annotated_prots else label
        node_labels[n] = "%s\n%d" % (label, ranks[n] if ranks[n] < first_zero_rank else first_zero_rank) if n in taxon_prots else label

        # maybe put the labels below the nodes?
        # helps with visualizing the background opacity
        graph_attr[n]['text-valign'] = 'bottom'
        # add the strain name to the popup
        attr_desc[n]['Strain'] = f_settings.TAX_TO_NAME[prot_species[n]]
        if n in predicted_prots:
            # don't need to normalize because the confidence values are already between 0 and 1
            if taxon and (n in non_taxon_annotated_prots or n in non_taxon_neg_prots):
                pass
            else:
                # UPDATE: use the node rank instead of the node score
                #graph_attr[n]['background-opacity'] = pred_local_conf[n]
                if n not in ranks:
                    graph_attr[n]['background-opacity'] = float(pred_local_conf[n])
                else:
                    #graph_attr[n]['background-opacity'] = float(pred_local_conf[n])
                    graph_attr[n]['background-opacity'] = max([0.9 - (ranks[n] / float(first_zero_rank)), float(pred_local_conf[n])])
                    attr_desc[n]["%s rank"%(alg_names[algorithm])] = ranks[n]
            if n in pred_cut_conf:
                attr_desc[n]["%s prediction score"%(algorithm)] = pred_cut_conf[n] 
            attr_desc[n]["%s prediction score"%(alg_names[algorithm])] = "%0.4f" % (float(pred_local_conf[n]))
        #elif n in annotated_prots or (taxon and (n in non_taxon_annotated_prots or n in non_taxon_neg_prots)) \
        #     or n in neg_prots:
            #if n in pred_local_conf:
            #    graph_attr[n]['background-opacity'] = pred_local_conf[n]
            #    attr_desc[n]["Local prediction confidence"] = pred_local_conf[n] 
        # also add the annotation to the popup
        if n in uniprot_to_evidencecode:
            codes = uniprot_to_evidencecode[n]
            # TODO add bullet points to the list
            #attr_desc[n]["Evidence code"] = ''.join(["%s (%s)\n" % (c, evidence_code_name[c]) for c in codes])
            # order it by exp, comp, then elec
            evidence_codes = ''.join(["<li>%s (%s)</li>" % (c, evidence_code_name[c]) for c in codes if evidence_code_type[c] == 'experimental'])
            evidence_codes += ''.join(["<li>%s (%s)</li>" % (c, evidence_code_name[c]) for c in codes if evidence_code_type[c] == 'computational'])
            evidence_codes += ''.join(["<li>%s (%s)</li>" % (c, evidence_code_name[c]) for c in codes if evidence_code_type[c] == 'electronic'])
            attr_desc[n]["Evidence code"] = "<ul>%s</ul>" % (evidence_codes)

    for u,v in prededges:
        e = (u,v)
        if e not in attr_desc:
            attr_desc[e] = {}
        if e not in graph_attr:
            graph_attr[e] = {}
        #attr_desc[e]["edge weight"] = G.adj[u][v]]['weight']
        if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version]:
            attr_desc[e]["Final edge weight"] = "%0.1f" % (W[node2idx[u]][:,node2idx[v]].A.flatten()[0])
            edge_type_weights = []
            # add the weights for the individual string networks
            for i in range(len(indices)):
                net_name = net_names[indices[i]]
                net_name = "SSN (E-value <= 0.1)" if 'seq-sim-e0_1' in net_name else net_name
                net = sparse_nets[indices[i]]
                w = net[node2idx[u]][:,node2idx[v]].A.flatten()[0]
                if w != 0:
                    #attr_desc[e][net_name] = "%0.1f" % (w)
                    edge_type_weights.append("<li>%s: %0.1f</li>" % (net_name, w))
            attr_desc[e]["Edge weights by type"] = "<ul>%s</ul>" % (''.join(sorted(edge_type_weights)))
        else:
            attr_desc[e]["Edge weight"] = "%0.1f" % (W[node2idx[u]][:,node2idx[v]].A.flatten()[0])
        # make the edges somewhat opaque for a better visual style
        graph_attr[e]['opacity'] = 0.7

    # set the width of the edges by the network weight
    edge_weights = {(u,v): float(W[node2idx[u]][:,node2idx[v]].A.flatten()[0]) for u,v in prededges}
    # TODO set the min and max as parameters or something
    #max_weight = 180 
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[opts.version]:
        max_weight = net_weights[0]*180
    else:
        max_weight = 180 
    for e in edge_weights:
        if edge_weights[e] > max_weight:
            edge_weights[e] = max_weight 
    graph_attr = gs.set_edge_width(prededges, edge_weights, graph_attr,
                                   a=1, b=12, min_weight=1, max_weight=max_weight)

    H = nx.Graph()
    H.add_edges_from(prededges)

    # see which DB the edge came from to set the edge color
    print("Getting the edge type from networks")
    if 'STRING' in f_settings.NETWORK_VERSION_INPUTS[version] and 'SEQ_SIM' in f_settings.NETWORK_VERSION_INPUTS[version]:
        print("\tFrom both STRING and SEQ_SIM")
        seq_sim_edges = set()
        for u,v in utils.readColumns(f_settings.SEQ_SIM_NETWORKS[version], 1, 2):
            #if (u,v) not in prededges:
            if not H.has_edge(u,v):
                continue
            # these are all undirected, so just store the sorted version
            u,v = tuple(sorted((u,v)))
            # give these the default color
            graph_attr[(u,v)]['color'] = edge_type_color['default']
            seq_sim_edges.add((u,v))

#        string_edges = set()
#        temp_version = '2017_10-string'
#        net = f_settings.NETWORK_template % (temp_version, temp_version) 
#        for u,v in utils.readColumns(net, 1, 2):
#            #if (u,v) not in prededges:
#            if not H.has_edge(u,v):
#                continue
#            # give these the default color
#            u,v = tuple(sorted((u,v)))
#            graph_attr[(u,v)]['color'] = edge_type_color['string']
#            string_edges.add((u,v))
        string_edges = prededges.difference(seq_sim_edges)
        print("\t%d edges from seq-sim, %d edges from STRING" % (len(seq_sim_edges), len(string_edges)))
        # set the color to STRING if it didn't come from sequence similarity
        for e in string_edges:
            #if 'color' not in graph_attr[e]:
            graph_attr[e]['color'] = edge_type_color['string']

    elif 'STRING' in f_settings.NETWORK_VERSION_INPUTS[version]:
        for e in graph_attr:
            graph_attr[e]['color'] = edge_type_color['string']
    else:
        for e in graph_attr:
            graph_attr[e]['color'] = edge_type_color['default']

    # apply the evidence code style to each protein
    for n in prednodes:
        if n in annotated_prots:
            graph_attr[n]['color'] = node_type_color['annotation']
        elif taxon and n in non_taxon_annotated_prots:
            graph_attr[n]['color'] = node_type_color['non-taxon-annotation']
        elif taxon and n in non_taxon_neg_prots:
            graph_attr[n]['color'] = node_type_color['non-taxon-neg-annotation']
        elif n in neg_prots:
            graph_attr[n]['color'] = node_type_color['neg-annotation']
        elif n in predicted_prots:
            graph_attr[n]['color'] = node_type_color['prediction']
        if n in uniprot_to_evidencecode:
            curr_style = ""
            for evidencecode in uniprot_to_evidencecode[n]:
                curr_type = evidence_code_type[evidencecode]
                if curr_type == "experimental":
                    curr_style = annotation_type_styles[curr_type]
                    break
                elif curr_style == "computational":
                    continue
                else:
                    curr_style = annotation_type_styles[curr_type]
            graph_attr[n].update(curr_style)

    # TODO build the popups here. That way the popup building logic can be separated from the
    # GSGraph building logic
    popups = {}
    prednodes = set([n for edge in prededges for n in edge])
    for n in prednodes:
        popups[n] = gs.buildNodePopup(n, attr_val=attr_desc)
    for u,v in prededges:
        popups[(u,v)] = gs.buildEdgePopup(u,v, node_labels=uniprot_to_gene, attr_val=attr_desc)

    # Now post to graphspace!
    print("Building GraphSpace graph")
    G = gs.constructGraph(prededges, node_labels=node_labels, graph_attr=graph_attr, popups=popups)

    # TODO add an option to build the 'graph information' tab legend/info
    # build the 'Graph Information' metadata
    #desc = gs.buildGraphDescription(opts.edges, opts.net)
    desc = ''
    metadata = {'description':desc,'tags':[], 'title':''}
    if tags is not None:
        metadata['tags'] = tags
    G.set_data(metadata)
    graph_name = "%s-%s-%s-%s-%s%s" % (selected_goname, selected_goid, algorithm, exp_name, version, postfix)
    G.set_name(graph_name)

    # rather than call it from here and repeat all the options, return G, and then call this after 
    #post_graph_to_graphspace(G, opts.username, opts.password, opts.graph_name, apply_layout=opts.apply_layout, layout_name=opts.layout_name,
    #                         group=opts.group, make_public=opts.make_public)
    return G, graph_name


# TODO modify the utils post_to_graphspace.py parseArgs so I can call it instead
def parseArgs(args):
    ## Parse command line args.
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)
    parser.add_option('','--version',type='string',
                      help="Version of the PPI to run. Options are: %s." % (', '.join(f_settings.ALLOWEDVERSIONS)))
    parser.add_option('-N','--net-file', type='string',
                     help="Network file to use. Default is the version's default network")
    parser.add_option('', '--exp-name', type='string',
                      help="Experiment name to use when running GAIN.")
    parser.add_option('', '--algorithm', type='string', metavar='STR',
                      help="Algorithm for which to get predictions. Options: '%s'" % ("', '".join(f_settings.ALGORITHM_OPTIONS.keys())))
    parser.add_option('', '--goid', type='string', metavar='STR',
                      help='GO-term ID for which annotations and precitions will be posted')
    parser.add_option('', '--pos-neg-file', type='string', action='append',
                      help="File containing positive and negative examples for each GO term")
    parser.add_option('-S', '--goid-summary-file', type='string', 
                      help="File containing GO term names and # of annotations for each GO term")
    parser.add_option('-T', '--taxon', type='string', 
                      help="Specify the species taxonomy ID to use. Otherwise, all species will be used")
    parser.add_option('', '--num-neighbors', type=int, default=1,
                      help="Number of neighbors to add around taxon positives. Default: 1")
    parser.add_option('', '--node-to-post', type='string', action="append",
                      help="UniProt ID of a taxon node for which to get neighbors. Can specify multiple")
    parser.add_option('', '--pos-neg-file-eval', type='string', action='append',
                     help="File containing positive and negative examples for each GO term used to evaluate for each species")
    parser.add_option('', '--keep-ann', action='store_true', default=False,
                     help="Don't leave out annotations when running the algorithms" +
                     "TODO allow for the option to run and evaluate everything together")
    parser.add_option('-a', '--alpha', type=float, default=1,
                     help="Alpha insulation parameter. Default=1")
    parser.add_option('', '--eps', type=float, default=0.0001,
                     help="Stopping criteria for SinkSource. Default: 0.0001")
    parser.add_option('', '--max-iters', 
                     help="Dummy option. Just here so it better matches the options of run_algs.py")
    # TODO add algorithm
    #parser.add_option('', '--edges', type='string', metavar='STR',
    #                  help='file edges to post.')
    #parser.add_option('', '--net', type='string', metavar='STR',
    #                  help='file of weighted directed edges. Can be used to get the weight of each edge for the popup. Tab-delimited file has columns TAIL,  HEAD,  WEIGHT.')
    #parser.add_option('', '--mapping-file', type='string', metavar='STR',
    #                  help='file used to map to a different namespace. Network/edge IDs (uniprot ids) should be in the first column with the other namespace (gene name) in the second')

    # posting options
    parser.add_option('-U', '--username', type='string', metavar='STR', 
                      help='GraphSpace account username to post graph to. Required')
    parser.add_option('-P', '--password', type='string', metavar='STR', 
                      help='Username\'s GraphSpace account password. Required')
    #parser.add_option('', '--graph-name', type='string', metavar='STR', default='test',
    #                  help='Graph name for posting to GraphSpace. Default = "test".')
    #parser.add_option('', '--outprefix', type='string', metavar='STR', default='test',
    #                  help='Prefix of name to place output files. Required.')
    parser.add_option('', '--name-postfix', type='string', metavar='STR', default='',
                      help='Postfix of graph name to post to graphspace.')
    parser.add_option('', '--group', type='string', metavar='STR', 
                      help='Name of group to share the graph with.')
    parser.add_option('', '--make-public', action="store_true", default=False,
                      help='Option to make the uploaded graph public')
    # TODO implement and test this option
    #parser.add_option('', '--group-id', type='string', metavar='STR',
    #                  help='ID of the group. Could be useful to share a graph with a group that is not owned by the person posting')
    parser.add_option('', '--tag', type='string', metavar='STR', action="append",
                      help='Tag to put on the graph. Can list multiple tags (for example --tag tag1 --tag tag2)')
    parser.add_option('', '--apply-layout', type='string', metavar='STR', 
                      help='Specify the name of a graph from which to apply a layout. Layout name specified by the --layout-name option. ' + 
                      'If left blank and the graph is being updated, it will attempt to apply the --layout-name layout.')
    parser.add_option('', '--layout-name', type='string', metavar='STR', default='layout1',
                      help="Name of the layout of the graph specified by the --apply-layout option to apply. Default: 'layout1'")
    # TODO implement parent nodes
    #parser.add_option('', '--include-parent-nodes', action="store_true", default=False,
    #                  help="Include source, target, intermediate parent nodes, and whatever parent nodes are in the --graph-attr file")

    # extra options
    #parser.add_option('', '--graph-attr', type='string', metavar='STR',
    #        help='Tab-separated File used to specify graph attributes 1st column: style, 2nd: style attribute, 3rd: list of uniprot ids (nodes/edges) separated by |, 4th: Description.')
    #parser.add_option('', '--set-edge-width', action="store_true", default=False,
    #                  help='Set edge widths according to the weight in the network file')

    (opts, args) = parser.parse_args(args)

    if opts.version not in f_settings.ALLOWEDVERSIONS:
        print("ERROR: '%s' not an allowed version. Options are: %s." % (opts.version, ', '.join(f_settings.ALLOWEDVERSIONS)))
        sys.exit(1)
    if opts.exp_name is None or opts.goid is None or opts.algorithm is None:
        print("--exp-name, --goid, --algorithm, required")
        sys.exit(1)

    if opts.algorithm not in set(f_settings.ALGORITHM_OPTIONS) \
       | {'localplus', 'local', 'sinksourceplus', 'sinksource', 'sinksource_squeeze', 'sinksourceplus_squeeze',
          'genemania', 'birgrank'}:
        print("--algorithm %s is not a valid algorithm name" % (opts.algorithm))
        sys.exit(1)

    if opts.username is None or opts.password is None:
        #parser.print_help()
        sys.exit("\nERROR: --username and --password required")

    if opts.net_file and not os.path.isfile(opts.net_file):
        sys.exit("\nERROR: --net-file %s does not exist" % (opts.net_file))

    return opts, args


if __name__ == '__main__':
    opts, args = parseArgs(sys.argv)
    # TODO allow for multiple versions. Will need to allow for multiple exp_names as well(?)
    version = opts.version
    G, graph_name = setup_post_to_graphspace(
        version, opts.algorithm, opts.exp_name,
        opts.goid, tags=opts.tag, pos_neg_files=opts.pos_neg_file,
        taxon=opts.taxon, goid_summary_file=opts.goid_summary_file,
        num_neighbors=opts.num_neighbors, nodes_to_post=opts.node_to_post,
        pos_neg_files_eval=opts.pos_neg_file_eval, keep_ann=opts.keep_ann,
        alpha=opts.alpha, eps=opts.eps, maxi=opts.max_iters,
    )
    if opts.taxon:
        graph_name += "-%s"%opts.taxon
    if opts.node_to_post is not None:
        graph_name += '-'.join(opts.node_to_post)
    graph_name += opts.name_postfix
    G.set_name(graph_name)
    # example command: python src/post_to_graphspace.py --edges toxic-sub-bind-edges.txt --net inputs/2017_10-seq-sim/2017_10-seq-sim-net.txt --username jeffl@vt.edu --password f1fan --graph-name toxic-substance-binding-cc-test5 --graph-attr graph_attr.txt --tag test1 --tag test2 --set-edge-width
    gs.post_graph_to_graphspace(G, opts.username, opts.password, graph_name, apply_layout=opts.apply_layout, layout_name=opts.layout_name,
                            group=opts.group, make_public=opts.make_public)
