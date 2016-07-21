"""
WARNING: Do not edit any files that this program relies on while it
is running. You will end up changing the experimental environment!

To safely run, copy this file and other files into another folder and 
run from that other location. Then you may edit the original files.

Usage: experiment.py <exp_folder>

Options:
  <exp_folder>          Path to destination folder for the sets of data.

"""

import docopt
from itertools import product
import os
from subprocess import check_call as cc
from tools import timeit

from analyze import summary
from batch import run_batch

def setup(batch_folder, methods, probs, flow_dict, force):

    def one_run(tup):
        species_depth, n_gene_trees, n_ind, Ne, n_sp, n_sp_trees = tup

        run_batch(batch_folder=batch_folder.format(*tup),
                  species_depth=species_depth,
                  n_gene_trees=n_gene_trees,
                  n_ind=n_ind,
                  n_sp=n_sp,
                  n_sp_trees=n_sp_trees,
                  Ne=Ne,
                  prob_missing=probs,
                  methods=methods,
                  flow_dict=flow_dict,
                  force=force)

    return one_run

if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    exp_folder = args['<exp_folder>']
    if exp_folder[0] != '/':
        exp_folder = os.path.abspath(os.path.join("..", exp_folder))

    if not os.path.isdir(exp_folder):
        os.makedirs(exp_folder)

    # c = [0.6, 1.0, 16.0]#0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    # c = list(map(float, c))
    # genes = [20]#, 40, 50, 60]
    # inds = [5]#, 6, 10]
    # methods = [1, 2]
    # probs = [0.1, 0.2, 0.05]
    # species = [3]#, 4, 5]
    # flow_dict = {"all":False,
    #              "generate":False,
    #              "drop":False,
    #              "impute":True,
    #              "analyze":False, 
    #              "--plus":True}

    # experiment run
    # c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 
    c = [1,2]
    c = list(map(float, c))
    genes = [10, 15]#, 40, 50, 60]
    inds = [2, 3]#, 6, 10]
    methods = [1, 2, 3, 4]
    probs = [0.1, 0.05]
    species = [2, 3]#, 4, 5]
    pop_size = [10000]
    depth = list(set(map(lambda x: int(x[0]*x[1]), product(c, pop_size))))
    trees = [2] # number of species trees
    flow_dict = {"all":True,
                 "generate":False,
                 "drop":False,
                 "impute":False,
                 "analyze":False, 
                 "--plus":True}
    force = False
    
    batch_folder = os.path.join(exp_folder, "d{}_g{}_i{}_n{}_s{}")

    f = lambda: list(map(setup(batch_folder, methods, probs, flow_dict, force),
                         product(depth, genes, inds, pop_size, species, trees)))
    timeit(f, "solving all problems")
    
    summary(exp_folder)
    cc("Rscript_$_summary.R_$_{}".format(exp_folder).split("_$_"))

