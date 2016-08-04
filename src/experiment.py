"""
Usage: experiment.py <exp_folder> [options]

Options:
  <exp_folder>          Path to destination folder for the sets of data.

  -f                    Overwrite files that may already exist within
                        the experiment folder.          [default: False]
  -t                    Testing mode.                   [default: False]
"""

import docopt
from itertools import product
import os
import subprocess

from analyze import compile_stats
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

    # get args
    exp_folder = args['<exp_folder>']
    test, experiment = args['-t'], not args['-t']
    force = args['-f']

    # set up filesystem
    if exp_folder[0] != '/':
        exp_folder = os.path.abspath(os.path.join("..", exp_folder))
    if not os.path.isdir(exp_folder):
        os.makedirs(exp_folder)

    # test run
    if test:
        c = [1]#, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
        c = list(map(float, c))
        genes = [10, 20]#, 40, 50, 60]
        inds = [5,10]#, 6, 10]
        methods = [1,2]#, 4]
        probs = [1]#, 0.05, 0.2]
        species = [3]#, 3, 5]
        pop_size = [10000]
        depth = list(set(map(lambda x: int(x[0]*x[1]), product(c, pop_size))))
        trees = [1] # number of species trees
        flow_dict = {"all": True,
                     "generate":False,
                     "drop":False,
                     "impute":False,
                     "analyze":True, 
                     "--plus":True}

    # experimental set up
    if experiment: 
        c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
        genes = [10, 20, 30]
        methods = [3, 4]
        probs = [0.05, 0.1, 0.2]
        species = [2, 3, 5]
        inds = [2, 5, 10]
        trees = [1]
        # c = [1, 4, 8, 12, 16, 20]
        # inds = [5, 10, 15]
        # methods = [1, 2]
        # probs = [4]
        # species = [3]
        # trees = [3] # number of species trees
        pop_size = [10000]
        c = list(map(float, c))
        depth = list(set(map(lambda x: int(x[0]*x[1]), product(c, pop_size))))
        
        # Options to set 
        flow_dict = {"all":False,        # overrides other options
                     "generate":False,  # generate trees
                     "drop":False,      # drop leaves
                     "impute":False,    # impute missing leaves
                     "analyze":True,   # analyze batches
                     "--plus":True}     # perform operations following
                                        # the first selection operation
    
    # batch folder naming scheme
    batch_folder = os.path.join(exp_folder, tools.batch_general)
    one_batch_run = setup(batch_folder, methods, probs, flow_dict, force)
    batch_iterator = product(depth, genes, inds, pop_size, species, trees)
    
    # choose run method depending on testing or not
    run_parallel= lambda: tools.parmap(one_batch_run, batch_iterator)
    run_serial = lambda: list(map(one_batch_run, batch_iterator))
    f = run_parallel if experiment else run_serial

    # run experiment
    tools.timeit(f, "solving all problems")
    
    # run statistics
    compile_stats(exp_folder)
    make_graphs = "Rscript_$_summary.R_$_{}".format(exp_folder)
    subprocess.check_call(make_graphs.split("_$_"))

