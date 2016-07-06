"""
Usage: generateTrees.py [options]

Options:
  -h, --help            Show this help message.
  -g, --n_gene_trees=.  Number of gene trees.            [default: 1000]
  -i, --n_ind=.         Number of individuals per species.  [default: 2]
  -n, --Ne=.            Effective population size.      [default: 10000]
  -o, --outfolder=.     Path to output folder.
                 [default: ../simulated_data/SD:Ne ratio sweep/py_nexus]
  -s, --n_sp=.          Number of species per species tree. [default: 5]
  -t, --n_sp_trees=.    Number of species trees.            [default: 2]

"""
import sys
import os
import string
import docopt
import dendropy as dp

# Some annoying but necessary abbreviations
taxa_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping
bd_tree = dp.simulate.birth_death_tree
cc_tree = dp.simulate.treesim.contained_coalescent_tree

def distance_from_root(node):
    """
    This function should exist but is not implemented cleanly in dendropy.
    """
    dist = 0
    while node.parent_node:
        dist += node.edge.length
        node = node.parent_node
    return dist

def species_tree(species, sp_depth):
    """
    Make a species tree and adjust the length of the edges to match the
    desired species depth.
    """
    t = bd_tree(birth_rate=1.0, death_rate=0, taxon_namespace=species)

    max_depth = max(map(distance_from_root, t.leaf_node_iter()))

    for e in t.preorder_edge_iter():
        e.length *= sp_depth/max_depth
    
    t.seed_node.edge.length = 0

    return t

def gen_trees(n_sp_trees, n_gene_trees, n_sp, n_ind, sp_depth, Ne=10000):
    # make taxa for species. names are "A", "B", "C", ...
    species = dp.TaxonNamespace(string.ascii_uppercase[:n_sp])

    # generate species trees and set population size of each edge to 10000
    # must explicitly make list, or cannot set pop_size
    sp_trees = list(map(lambda x: species_tree(species, sp_depth),
                        range(n_sp_trees)))
    for tree in sp_trees:
        for edge in tree.postorder_edge_iter():
            setattr(edge, 'pop_size', Ne)

    # convert species names to individual names and build map between taxa
    label = lambda taxon, index: "{}{}".format(taxon.label.lower(), index + 1)
    si_map = taxa_map(containing_taxon_namespace=species,
                      num_contained=n_ind,
                      contained_taxon_label_fn=label)

    # make contained coalescent trees
    make_ctrees = lambda tree: dp.TreeList(map(lambda y: cc_tree(tree, si_map),
                                               range(n_gene_trees)))
    gene_trees = list(map(make_ctrees, sp_trees))

    return dp.TreeList(sp_trees), gene_trees

if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    print(args)
    n_sp_trees = int(args['--n_sp_trees'])
    n_gene_trees = int(args['--n_gene_trees'])
    n_sp = int(args['--n_sp'])
    n_ind = int(args['--n_ind'])
    Ne = int(args['--Ne'])
    outfolder = args['--outfolder']

    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    gene_formula = 'c{}_genes_{}.nex'
    sp_formula = 'c{}_species.nex'

    for c in [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]:
        sp_depth = Ne*c

        sp, gn = gen_trees(n_sp_trees, n_gene_trees, n_sp, n_ind, sp_depth)

        sp_path = os.path.join(outfolder, sp_formula.format(c))
        sp.write(path=sp_path, schema="nexus")
        
        for i, treelist in enumerate(gn):

            outfile = os.path.join(outfolder, gene_formula.format(c, i+1))
            treelist.write(path=outfile, schema="nexus")

