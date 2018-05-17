# TODO The deepnf code is at /home/jeffl/git-workspace/from-lit/deepNF

print("Importing libraries")
import os, sys
#import matlab.engine
#from optparse import OptionParser
# add the folder above to the path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils.file_utils as utils
# add two folders up for fungcat settings
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import fungcat_settings as f_settings
import networkx as nx
import pandas as pd
import numpy as np
from scipy.io import savemat
from scipy import sparse
from tqdm import tqdm
#import igacat.go_term_prediction_examples.go_term_prediction_examples as go_examples


def convert_nodes_to_int(G):
    index = 0
    node2int = {}
    int2node = {}
    for n in G.nodes():
        node2int[n] = index
        int2node[index] = n
        index += 1
    # see also convert_node_labels_to_integers
    G = nx.relabel_nodes(G,node2int, copy=False)
    return G, node2int, int2node


out_dir = "test/deepnf"
utils.checkDir(out_dir)

#version = "2017_10-seq-sim-x5-string"
version = "2017_10-string"
# start with a single species
taxon = "208964"

#prots_file = "%s/%s-prots.txt" % (out_dir, version)
#out_file = "%s/%s-net.mat" % (out_dir, version)
prots_file = "%s/%s-%s-prots.txt" % (out_dir, version, taxon)
out_file = "%s/%s-%s-adj-net.txt" % (out_dir, version, taxon)
if os.path.isfile(prots_file) and os.path.isfile(out_file):
    print("Skipping mapping nodes to ids. %s exists" % (out_file,))
    #prots = utils.readItemList(prots_file)
else:
    # G = readtable(filename)
    network_file = "inputs/2017_10-seq-sim-x5-string/2017_10-seq-sim-x5-string-net.txt"
    if taxon is not None:
        network_file = f_settings.STRING_TAXON_UNIPROT % (taxon, taxon, f_settings.STRING_CUTOFF)
    print("Reading network from %s" % (network_file))


    G = nx.Graph()
    #G.add_weighted_edges_from([(u,v,float(w)) for u,v,w in lines])
    #edges = {}
    #nodes = set()
    with open(network_file, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            u,v,w = line.rstrip().split('\t')[:3]
            G.add_edge(u,v,weight=float(w))
            #edges[(u,v)] = float(w)
            #nodes.update(set([u,v]))
    print("\t%d nodes and %d edges" % (G.number_of_nodes(), G.number_of_edges()))

    print("\trelabeling node IDs with integers")
    G, node2int, int2node = convert_nodes_to_int(G)
    print("\twriting network to %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write(''.join(["%s\t%s\t%s\n" % (u,v,data['weight']) for u,v,data in G.edges(data=True)]))
    print("\twriting node2int labels to %s" % (prots_file))
    with open(prots_file, 'w') as out:
        out.write(''.join(["%s\t%s\n" % (int2node[n], n) for n in sorted(int2node)]))

print("Setting up the annotation matrix not yet implemented. Quitting.")
sys.exit()

# now read the annotations and map those to IDs
#ann_file = "%s/%s-%s-%s-goa.txt" % (out_Dir, version, taxon)
#if os.path.isfile(ann_file):
#    print("skipping"
#else:
#print("Mapping go annotations to node IDs to %s" % (ann_file))
#positives = np.asarray([node2int[n] for n in all_goid_prots[goterm] if n in node2int])



# this file just has the direct GO annotations (not propagated)
#gaf_file = "/data/inputs/goa/taxon/19-strains-goa.gaf"
# this file has the propagated annotations for all 19 species
all_taxon_goa = '/data/inputs/goa/taxon/all-taxon-goa.txt'
gain_file = all_taxon_goa
if taxon is not None:
    gain_file = f_settings.FUN_FILE % (taxon, taxon)

# columns: 'orf', 'goid', 'hierarchy', 'evidencecode', 'annotation type' (described here: http://bioinformatics.cs.vt.edu/~murali/software/biorithm/gain.html)
print("reading annotations from %s" % (gain_file))
df = pd.read_csv(gain_file, sep='\t')
print("%d proteins have at least 1 annotation" % (len(set(df['orf'].tolist()))))
print("%d goids have at least 1 protein annotated to it" % (len(set(df['goid'].tolist()))))
df['hierarchy'] = df['hierarchy'].apply(lambda x: x.upper())
# also add the leading GO:000 to the ID
df['goid'] = df['goid'].apply(lambda x: "GO:" + "0"*(7-len(str(x))) + str(x))

#for h in ["C", "F", "P"]:
for h in ["P"]:
    #print(h)
    df_h = df[df["hierarchy"] == h]
    df_h = df_h[["orf", "goid"]]
    #prots = sorted(set(df_h['orf'].tolist()))
    goids = sorted(set(df_h['goid'].tolist()))
    print("%d goids have at least 1 %s annotation" % (len(goids), h))

    # write this order of the GO IDs to a file so they can easily be accessed (to get the rows of the sparse matrix)
    # without having to perform these operations again
    out_file = "%s/%s-goids.txt" % (out_dir, h)
    print("Writing %s" % (out_file))
    with open(out_file, 'w') as out:
        out.write('\n'.join(goids)+'\n')

    # now convert it to a table and write it to a file
    # convert a categorical variable into a "dummy" or "indicator" DataFrame
    # moves the column to the header and puts a 0 or 1 for if the value is present or not
    # UPDATE: this hits a memory limit. build the table manually instead
    #annotation_matrix = pd.get_dummies(df_h['goid']).groupby(df_h['orf']).apply(max)

    prot_goids = {orf: set(goids['goid'].tolist()) for orf, goids in df_h[['orf', 'goid']].groupby('orf')}
    # TODO figure out a better way to build this matrix
    # UPDATE: don't write the header line and row names for the sparse matrix
    #table = [[''] + goids]
    print("Building an annotation table where rows are GO terms and columns are proteins")
    table = []
    for prot in tqdm(prots):
        #row = [prot]
        row = []
        curr_goids = prot_goids.pop(prot, set())
        for goid in goids:
            if goid in curr_goids:
                row.append(1)
            else:
                row.append(0)
        table.append(row)

    # convert it to a sparse matrix and transpose it so it matches the GO DAG and G network (above)
    print("Converting the annotation table to a sparse matrix")
    annotation_matrix = sparse.csr_matrix(table).transpose()
    #annotation_matrix.to_csv(out_file)
    #with open(out_file, 'w') as out:
    #    for row in table:
    #        out.write(','.join(row) + '\n')

    # convert the GO DAG to a sparse matrix, while maintaining the order of goids so it matches with the annotation matrix
    #dag_matrix = nx.to_scipy_sparse_matrix(G, nodelist=goids, weight=None)
    out_file = "%s/%s-annotations.mat" % (out_dir, h)
    print("writing %s" % (out_file))
    savemat(out_file, annotation_matrix)
    #savemat(out_file, {'R':annotation_matrix, 'H': dag_matrix})
