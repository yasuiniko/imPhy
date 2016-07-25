"""
Usage: experiment.py <exp_folder>

Options:
  <exp_folder>          Path to destination folder for the sets of data.

"""

import docopt
from itertools import product
import os
from subprocess import check_call as cc

from analyze import summary
from batch import run_batch
import tools

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

    # # test run
    # c = [0.6]#, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    # c = list(map(float, c))
    # genes = [10, 20, 30]#, 40, 50, 60]
    # inds = [2, 5, 10]#, 6, 10]
    # methods = [3]#, 4]
    # probs = [0.1]#, 0.05, 0.2]
    # species = [2]#, 3, 5]
    # pop_size = [10000]
    # depth = list(set(map(lambda x: int(x[0]*x[1]), product(c, pop_size))))
    # trees = [1] # number of species trees
    # flow_dict = {"all":True,
    #              "generate":False,
    #              "drop":False,
    #              "impute":False,
    #              "analyze":False, 
    #              "--plus":True}
    # force = True

    # experiment run
    c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    c = list(map(float, c))
    genes = [10, 20, 30]#, 40, 50, 60]
    inds = [2, 5, 10]#, 6, 10]
    methods = [3, 4]
    probs = [0.1, 0.05, 0.2]
    species = [2, 3, 5]
    pop_size = [10000]
    depth = list(set(map(lambda x: int(x[0]*x[1]), product(c, pop_size))))
    trees = [1] # number of species trees
    flow_dict = {"all":False,
                 "generate":False,
                 "drop":False,
                 "impute":False,
                 "analyze":True, 
                 "--plus":True}
    force = True
    
    batch_folder = os.path.join(exp_folder, "d{}_g{}_i{}_n{}_s{}")
    one_batch_run = setup(batch_folder, methods, probs, flow_dict, force)
    batch_iterator = product(depth, genes, inds, pop_size, species, trees)
    f = lambda: tools.parmap(one_batch_run, batch_iterator)
    tools.timeit(f, "solving all problems")
    
    summary(exp_folder)
    cc("Rscript_$_summary.R_$_{}".format(exp_folder).split("_$_"))

