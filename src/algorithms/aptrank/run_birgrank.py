
import sys, os
from optparse import OptionParser, OptionGroup
sys.path.append('src')
import fungcat_settings as f_settings
import utils.file_utils as utils
import algorithms.setup_sparse_networks as ssn
import algorithms.gain_scipy.alg_utils as alg_utils
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from birgrank import birgRank
import igacat.go_term_prediction_examples.go_term_prediction_examples as go_examples
from scipy import sparse
from scipy.io import savemat, loadmat
from networkx.algorithms.dag import descendants
import networkx as nx


def parse_args(args):
    usage = '%s [options]\n' % (sys.argv[0])
    parser = OptionParser(usage=usage)

    group = OptionGroup(parser, 'Main Options')
    group.add_option('-N','--net-file',type='string',
                     help="Network file to use.")
    group.add_option('', '--pos-neg-file', type='string', action='append',
                     help="File containing positive and negative examples for each GO term. Must use either this or the --gaf-file option")
    group.add_option('-g', '--gaf-file', type='string',
                     help="File containing GO annotations in GAF format. Annotations will not be propagated")
    group.add_option('-b', '--obo-file', type='string', default=f_settings.GO_FILE,
                     help="GO OBO file which contains the GO DAG. Default: %s" % (f_settings.GO_FILE))
    group.add_option('-H', '--hierarchy', type='string', default='bp',
                     help="Hierarchy to use when creating the matrices and running BirgRank. Default: 'bp'")
    group.add_option('', '--ignore-ec', type='string',
                     help="Comma-separated list of evidence codes where annotations with the specified codes will be ignored when parsing the GAF file. " +
                     "For example, specifying 'IEA' will skip all annotations with an evidence code 'IEA'. ")
    group.add_option('-o', '--out-pref', type='string',
                     help="Output prefix for the annotations matrix and heirarchy matrix to create. " +
                     "Default: %s" % ('TODO'))
    parser.add_option_group(group)

    #group = OptionGroup(parser, 'BirgRank Options')
    #group.add_option('-a', '--alpha', type=float, action="append",
    #                 help="Alpha insulation parameter. Default=0.8")
    #group.add_option('-W', '--num-pred-to-write', type='int', default=100,
    #                 help="Number of predictions to write to the file. If 0, none will be written. If -1, all will be written. Default=100")
    #group.add_option('', '--only-cv', action="store_true", default=False,
    #                 help="Perform cross-validation only")
    #group.add_option('-C', '--cross-validation-folds', type='int',
    #                 help="Perform cross validation using the specified # of folds. Usually 5")

    (opts, args) = parser.parse_args(args)
    return opts


# take in a network, and annotation matrix
# also take in a GO heirarchy and make it into a matrix
def main(sparse_net_file, obo_file, pos_neg_files=None, gaf_file=None, ignore_ec=["IEA"],
         alpha=.5, theta=.5, mu=.5, h="bp", out_pref=None):
#):
    print("Reading network from %s" % (sparse_net_file))
    sparse_net = sparse.load_npz(sparse_net_file)
    node2idx_file = sparse_net_file + "-node-ids.txt"
    print("Reading node names from %s" % (node2idx_file))
    if not os.path.isfile(node2idx_file):
        node2idx_file = sparse_net_file.replace('-normalized.npz', '-node2idx.txt')
        print("\tDoesn't exist. Checking %s" % (node2idx_file))
    prots = utils.readItemList(node2idx_file, 1)

    # parse the go_dags first as it also sets up the goid_to_category dictionary
    go_dags = go_examples.parse_obo_file_and_build_dags(obo_file)

    dag_matrix, an_matrix, goids = build_h_ann_matrices(prots, go_dags, pos_neg_files=None, gaf_file=None, h='bp')
    # make sure they're type float so matlab will parse them correctly
    sparse_net = sparse_net.astype('float') 
    ann_matrix = ann_matrix.astype('float') 
    dag_matrix = dag_matrix.astype('float')

    if out_pref is not None:
        out_file = "%s%s-annotations-and-go-dag.mat" % (out_pref, h)
        utils.checkDir(os.path.dirname(out_file))

        print("\twriting graph, annotation, and hierarchy matrices to %s" % (out_file))
        # write these to a file to run the matlab BirgRank 
        savemat(out_file, {"G": sparse_net, "R": ann_matrix, "H": dag_matrix}, do_compression=True)

        goids_file = "%s%s-goids.txt" % (out_pref, h)
        print("\twriting goids to %s" % (goids_file))
        with open(goids_file, 'w') as out:
            out.write(''.join("%s\n" % (goid) for goid in goids))

    run_birgrank = True 
    if run_birgrank is True:
        Xh = birgRank(sparse_net, ann_matrix.transpose(), dag_matrix, alpha=.5, theta=.5, mu=.5)

        # TODO write the results to a file (?)
        # for now, just print a single value to see if it matches with Matlab
        print(Xh[55986, 1])
    return


def build_h_ann_matrices(
        prots, go_dags, pos_neg_files=None, gaf_file=None, h='bp',
        goterms=None):
    """
    """
    category_to_name = {"C": "cc", "P": "bp", "F": "mf"}
    name_to_category = {val: key for key, val in category_to_name.items()}
    cat = name_to_category[h]

    # mapping from a GO term ID and the category it belongs to ('C', 'F' or 'P')
    goid_to_category = go_examples.goid_to_category

    # aptrank doesn't use negatives, so just get the positives
    if pos_neg_files is not None:
        goid_prots, _ = alg_utils.parse_pos_neg_files(pos_neg_files, goterms=goterms) 
    elif gaf_file is not None:
        prot_goids_by_c, goid_prots, _, _ = go_examples.parse_gaf_file(
            gaf_file, ignore_ec=ignore_ec) 
        goid_prots = {goid: prots for goid, prots in goid_prots.items()
                      if goid in goid_to_category and goid_to_category[goid] == cat}
    else:
        print("ERROR: must specify a pos-neg-file or gaf-file")
        sys.exit(1)
    print("\t%d GO terms" % (len(goid_prots)))
    if len(goid_prots) == 0:
        print("\tskipping")
        return sparse.csr_matrix((0,0)), sparse.csr_matrix((0,0)), []

    dag_matrix, goids = build_hierarchy_matrix(go_dags[cat], goid_prots.keys(), h=h)
    ann_matrix = setup_sparse_annotations(goid_prots, prots, goids, h=h)

    return dag_matrix, ann_matrix, goids


# TODO could also just use the annotation matrix and keep the entries > 0
# however, because AptRank uses the hierarchy,
# it should keep all of the GO terms higher up and not be cutoff (e.g., those with < 1000)
# However, I should probably also only use the leaf terms from the GAF file
def setup_sparse_annotations(
        goid_prots, prots, goids, h="bp", taxon=None):
    """
    *goid_prots*: dictionary containing the proteins annotated to each go term
    *prots*: list of proteins in the same order as the matrix
    *goids*: list of goids in the same order as the hierarchy matrix
    *returns*: A matrix with goterm rows, protein/node columns, and 1,0 for pos,unk values
    """
    node2idx = {prot: i for i, prot in enumerate(prots)}
    goid2idx = {goid: i for i, goid in enumerate(goids)}
    # rather than explicity building the matrix, use the indices to build a coordinate matrix
    # rows are prots, cols are goids
    i_list = []
    j_list = []
    data = []
    # limit the annotations to the proteins which are in the networks
    prots_set = set(prots)
    print("%d goids in Hierarchy, %d goids with annotations" %
          (len(goids), len(goid_prots)))
    for goid in goid_prots:
        for prot in goid_prots[goid] & prots_set:
            i_list.append(node2idx[prot])
            j_list.append(goid2idx[goid])
            data.append(1)

#    # limit it to the current taxon
#    # used if the prots are from multiple species (as is the case in the pos_neg_files
#    if taxon is not None:
#        print("Getting species of each prot from %s" % (f_settings.UNIPROT_TO_SPECIES))
#        #print("Limiting the prots to those for taxon %s (%s)" % (taxon, selected_species[taxon]))
#        print("Limiting the prots to those for taxon %s" % (taxon))
#        uniprot_to_species = utils.readDict(f_settings.UNIPROT_TO_SPECIES, 1,2)
#        # also build the reverse
#        species_to_uniprot = defaultdict(set)
#        for p in uniprot_to_species:
#            species_to_uniprot[uniprot_to_species[p]].add(p)
#        if taxon not in species_to_uniprot:
#            print("Error: taxon ID '%d' not found" % (taxon))
#            sys.exit()
#        taxon_prots = species_to_uniprot[taxon]
#        # also limit the proteins to those in the network
#        print("\t%d prots for taxon %s. Limiting to the %d in the network" % (len(taxon_prots), taxon, len(node2idx)))
#        taxon_prots = taxon_prots & set(node2idx.keys())
#        goid_prots = {goid: prots & taxon_prots for goid, prots in goid_prots.items()}
#        goids_to_remove = set()
#        for goid in goid_prots:
#            if len(goid_prots[goid]) == 0:
#                goids_to_remove.add(goid)
#        if len(goids_to_remove) > 0:
#            print("Warning: %d goterms have 0 annotations. Removing them." % (len(goids_to_remove)))
#            for goid in goids_to_remove:
#                goid_prots.pop(goid)

    # convert it to a sparse matrix 
    print("Building a sparse matrix of annotations")
    ann_matrix = sparse.coo_matrix((data, (i_list, j_list)), shape=(len(prots), len(goids)), dtype=float).tocsr()
    print("\t%d pos/neg annotations" % (len(ann_matrix.data)))
    ann_matrix = ann_matrix.transpose()

    return ann_matrix


def build_hierarchy_matrix(go_dag, goids, h="bp"):
    """
    *goids*: the leaf terms to use to get a sub-graph of the DAG.
        All ancestor terms will be included in the DAG
    """

    # UPDATE: limit to only the GO terms in R
    print("Limiting DAG to only the %d %s GO terms that have at least 1 annotation (assuming annotations already propagated up the DAG)" % (len(goids), h))
    ancestor_goids = set()
    for goid in goids:
        # if we already have the ancestors of this goid, then skip
        if goid in ancestor_goids:
            continue
        ancestor_goids.update(descendants(go_dag, goid))
    ancestor_goids.update(goids)
    goids_list = sorted(ancestor_goids)

    G = nx.subgraph(go_dag, ancestor_goids)
    print("\t%s DAG has %d nodes and %d edges" % (h, G.number_of_nodes(), G.number_of_edges()))

    # convert the GO DAG to a sparse matrix, while maintaining the order of goids so it matches with the annotation matrix
    dag_matrix = nx.to_scipy_sparse_matrix(G, nodelist=goids_list, weight=None)
    #out_file = "%s/%s-annotations-and-go-dag.mat" % (out_dir, h)
    #print("writing %s" % (out_file))
    #savemat(out_file, {'R':annotation_matrix, 'H': dag_matrix})

    return dag_matrix, goids_list


if __name__ == "__main__":
    opts = parse_args(sys.argv)
    ignore_ec = [] if opts.ignore_ec is None else opts.ignore_ec.split(',') 
    main(opts.net_file, opts.obo_file, opts.pos_neg_file, opts.gaf_file,
         ignore_ec=ignore_ec, out_pref=opts.out_pref)
