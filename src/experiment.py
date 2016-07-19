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

def setup(batch_folder, methods, probs, flow_dict):

    def one_run(tup):
        c, n_gene_trees, n_ind, n_sp = tup

        run_batch(batch_folder=batch_folder.format(*tup),
                  c=c,
                  n_gene_trees=n_gene_trees,
                  n_ind=n_ind,
                  n_sp=n_sp,
                  prob_missing=probs,
                  Ne=10000,
                  methods=methods,
                  flow_dict=flow_dict)

    return one_run

if __name__ == "__main__":
    args = docopt.docopt(__doc__)
    exp_folder = args['<exp_folder>']
    if exp_folder[0] != '/':
        exp_folder = os.path.abspath(os.path.join("..", exp_folder))

    if not os.path.isdir(exp_folder):
        os.makedirs(exp_folder)

    # edit these lists to your heart's desire
    # c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    # c.reverse()
    # genes = [10, 20, 30, 40, 50, 60]
    # inds = [4, 5, 6, 10]
    # methods = [1, 2]
    # probs = [0.1, 0.2, 0.5]
    # species = [2, 3, 4, 5]


    c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    c = list(map(float, c))
    genes = [10, 20, 30]#, 40, 50, 60]
    inds = [2, 5, 10]#, 6, 10]
    methods = [1, 2]
    probs = [0.1, 0.2, 0.05]
    species = [2, 3, 5]#, 4, 5]
    flow_dict = {"all":True,
                 "generate":False,
                 "drop":False,
                 "impute":False,
                 "analyze":False, 
                 "--plus":True}

    # # testing
    # c = [0.6]
    # genes = [10]
    # inds = [4]
    # methods = [1]
    # probs = [0.1]
    # species = [2]

    # names of the fields used in the format string, in the same order
    # as they appear in the itertools.product tuple.
    names = ["c", "genes", "inds", "sp"]
    
    # feel free also to change the experiment batch path
    batch_folder = os.path.join(exp_folder, "c{}_g{}_i{}_s{}")

    f = lambda: list(map(setup(batch_folder, methods, probs, flow_dict),
                         product(c, genes, inds, species)))
    timeit(f, "solving all problems")
    summary(exp_folder)
    cc("Rscript_$_summary.R_$_{}".format(exp_folder).split("_$_"))

