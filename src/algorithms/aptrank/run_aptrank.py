
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
from scipy.io import savemat
from scipy import sparse
from tqdm import tqdm
import igacat.go_term_prediction_examples.go_term_prediction_examples as go_examples

out_dir = "test/aptrank2"
utils.checkDir(out_dir)

#version = "2017_10-seq-sim-x5-string"
version = "2017_10-string"
# start with a single species
taxon = "208964"

#prots_file = "%s/%s-prots.txt" % (out_dir, version)
#out_file = "%s/%s-net.mat" % (out_dir, version)
prots_file = "%s/%s-%s-prots.txt" % (out_dir, version, taxon)
out_file = "%s/%s-%s-net.mat" % (out_dir, version, taxon)
if os.path.isfile(prots_file) and os.path.isfile(out_file):
    print("Skipping generating the network matrix %s. Reading protein ordering from %s" % (out_file, prots_file))
    prots = utils.readItemList(prots_file)
else:
    # G = readtable(filename)
    network_file = "inputs/2017_10-seq-sim-x5-string/2017_10-seq-sim-x5-string-net.txt"
    if taxon is not None:
        network_file = f_settings.STRING_TAXON_UNIPROT % (taxon, taxon, f_settings.STRING_CUTOFF)
    print("Reading network from %s" % (network_file))
    lines = utils.readColumns(network_file, 1,2,3)
    G = nx.Graph()
    # convert the weights to be betwee 0 and 1
    G.add_weighted_edges_from([(u,v, float(w)/1000.0) for u,v,w in lines])
    #print(G.edges['Q9I391','Q9HZJ2'])
    # write the network as an adjacency matrix table
    prots = G.nodes()

    # write the prots to a file
    print("Writing %s" % (prots_file))
    with open(prots_file, 'w') as out:
        out.write('\n'.join(prots)+'\n')

    print("Converting graph to an adjaceny table")
    # now build the adjacency matrix as a table
    print("Writing the sequence similarity + string network as an adjacency matrix to %s" % (out_file))
    #with open(out_file, 'w') as out:
    #    out.write(','.join([''] + prots) + '\n')
    #    for prot in tqdm(prots):
    #        row = [prot]
    #        for prot2 in prots:
    #            if G.has_edge(prot, prot2):
    #                row.append(G.edges[prot, prot2]['weight'])
    #            else:
    #                row.append('0')
    #        out.write(','.join(row) + '\n')

    # default nodelist order is G.nodes()
    sparse_matrix = nx.to_scipy_sparse_matrix(G, nodelist=prots)
    # save this graph into its own matlab file because it doesn't change by GO category
    savemat(out_file, {'G':sparse_matrix})

#opts, args = go_term_prediction_examples.parse_args(sys.argv) 
# build the GO DAG
# keeps only the is_a relationships for now
obo_file = "/data/inputs/goa/2017-09-26-go.obo"
go_dags = go_examples.parse_obo_file_and_build_dags(obo_file)

# write a matrix out of the DAG
# just use the biological process for now
#for c in go_dags:
    #nx.write_adjlist(go_dags[c], out_file)
    #sparse_matrix = nx.adjacency_matrix(go_dags[c])
    #savemat(out_file, {'H':sparse_matrix})
    # this is way too huge
    #df = nx.to_pandas_dataframe(go_dags[c]).astype(int)
    #out_file = "test/%s-go-dag.csv" % (c)
    #print("writing %s" % (out_file))
    #df.to_csv(out_file)

    # just write it to an edge list and make it into a matrix using matlab
    # UPDATE: limit to only the GO terms in R
#    out_file = "test/%s-go-dag.tsv" % (c)
#    print("writing %s" % (out_file))
#    nx.write_edgelist(go_dags[c], out_file, delimiter='\t', data=False)

# now build a matrix out of the annotations
# 0 or 1 for if an annotation is present or not
#gaf_file = "/data/inputs/goa/taxon/19-strains-goa.gaf"
#gain_file = "/data/inputs/goa/taxon/208964/208964-goa-all-fun.txt"
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

    # UPDATE: limit to only the GO terms in R
    print("Limiting DAG to only the %d %s GO terms that have at least 1 annotation (assuming annotations already propagated up the DAG)" % (len(goids), h))
    G = nx.subgraph(go_dags[h], goids)
    print("\t%s DAG has %d nodes and %d edges" % (h, G.number_of_nodes(), G.number_of_edges()))

    # convert the GO DAG to a sparse matrix, while maintaining the order of goids so it matches with the annotation matrix
    dag_matrix = nx.to_scipy_sparse_matrix(G, nodelist=goids, weight=None)
    out_file = "%s/%s-annotations-and-go-dag.mat" % (out_dir, h)
    print("writing %s" % (out_file))
    savemat(out_file, {'R':annotation_matrix, 'H': dag_matrix})
#    # now build the adjacency matrix as a table
#    table = [[''] + goids]
#    for goid in tqdm(goids):
#        row = [goid]
#        for goid2 in goids:
#            if G.has_edge(goid, goid2):
#                row.append('1')
#            else:
#                row.append('0')
#        table.append(row)
#    out_file = "%s/%s-go-dag.csv" % (out_dir, h)
#    print("Writing hierarchy as an adjacency matrix to %s" % (out_file))
#    with open(out_file, 'w') as out:
#        for row in table:
#            out.write(','.join(row) + '\n')
