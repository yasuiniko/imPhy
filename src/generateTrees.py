"""
Usage: generateTrees.py [<d>... options]

Options:
  -h, --help            Show this help message.
  <d>                   Species depths to consider
  -g, --n_gene_trees=.  Number of gene trees.            [default: 1000]
  -i, --n_ind=.         Number of individuals per species.  [default: 2]
  -n, --Ne=.            Effective population size.      [default: 10000]
  -o, --outfolder=.     Path to output folder.
                 [default: ../simulated_data/SD:Ne ratio sweep/py_nexus]
  -s, --n_sp=.          Number of species.                  [default: 5]
  -t, --n_sp_trees=.    Number of species trees.            [default: 1]

"""

import decimal
import dendropy as dp
import docopt
import os
import string
import sys

import tools

# Some annoying but necessary abbreviations
taxa_map = dp.TaxonNamespaceMapping.create_contained_taxon_mapping
bd_tree = dp.simulate.birth_death_tree
cc_tree = dp.simulate.treesim.contained_coalescent_tree

def distance_from_root(node):
    """
    This function should exist but is not implemented cleanly in dendropy.
    """
    dist = 0 #Decimal('0')
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
    
    # # convert distances to Decimals
    # for e in t.preorder_edge_iter():
    #         e.length = Decimal(str(e.length))
    
    max_depth = max(map(distance_from_root, t.leaf_node_iter()))

    # scale distances
    for e in t.preorder_edge_iter():
        # e.length *= Decimal(str(sp_depth))/max_depth
        e.length *= sp_depth/max_depth

    # # convert back to float
    # for e in t.preorder_edge_iter():
    #     e.length = float(e.length)
    
    t.seed_node.edge.length = 0

    return t

def gen_trees(n_sp_trees, n_gene_trees, n_sp, n_ind, sp_depth, Ne):
    # make taxa for species. names are "A", "B", "C", ...
    species = dp.TaxonNamespace(string.ascii_uppercase[:n_sp])

    # generate species trees and set population size of each edge to 10000
    # must explicitly make list, or cannot set pop_size
    sp_trees = dp.TreeList(map(lambda x: species_tree(species, sp_depth),
                               range(n_sp_trees)),
                           taxon_namespace=species)
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

    return sp_trees, gene_trees

def generateTrees(d_list,
                  n_gene_trees,
                  n_sp,
                  n_ind,
                  Ne,
                  n_sp_trees,
                  outfolder):
    # set up filesystem
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)
    gene_formula = tools.gene
    sp_formula = 'species.nex'

    for sp_depth in d_list:

        sp, gn = gen_trees(n_sp_trees, n_gene_trees, n_sp, n_ind, sp_depth, Ne)

        sp_path = os.path.join(outfolder, sp_formula)
        sp.write(path=sp_path, schema="nexus")
        
        for i, treelist in enumerate(gn):

            outfile = os.path.join(outfolder, tools.gene.format(sp_depth,
                                                                n_gene_trees,
                                                                n_ind,
                                                                Ne,
                                                                n_sp,
                                                                i+1))
            treelist.write(path=outfile, schema="nexus")

if __name__ == "__main__":

    # get args
    args = docopt.docopt(__doc__)
    n_sp_trees = int(args['--n_sp_trees'])
    n_gene_trees = int(args['--n_gene_trees'])
    n_sp = int(args['--n_sp'])
    n_ind = int(args['--n_ind'])
    Ne = int(args['--Ne'])
    outfolder = args['--outfolder']

    generateTrees(args['<d>'],
                  n_gene_trees,
                  n_sp,
                  n_ind,
                  Ne,
                  n_sp_trees,
                  outfolder)



