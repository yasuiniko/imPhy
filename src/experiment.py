"""
Usage: experiment.py <exp_folder> [options]

Options:
  <exp_folder>          Path to destination folder for the sets of data.

  -f                    Overwrite files that may already exist within
                        the experiment folder.          [default: False]
  -p                    Parallel mode.                  [default: False]
  -t                    Testing mode.                   [default: False]
"""

import docopt
from functools import partial
from itertools import product
import os
import shutil
import subprocess

from compile_stats import compile_stats
from batch import run_batch
import tools

def setup(batch_folder, methods, probs, flow_dict, force, tup):

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

if __name__ == "__main__":
    args = docopt.docopt(__doc__)

    # get args
    exp_folder = args['<exp_folder>']
    test, experiment = args['-t'], not args['-t']
    parallel = args['-p']
    force = args['-f']

    # set up filesystem
    if exp_folder[0] != '/':
        exp_folder = os.path.abspath(os.path.join("..", exp_folder))
    if not os.path.isdir(exp_folder):
        os.makedirs(exp_folder)
    heatpath = os.path.join(exp_folder, 'heatmaps')
    if os.path.isdir(heatpath):
        shutil.rmtree(heatpath)
    os.makedirs(heatpath)

    # test run
    if test:
        c = [20]#, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
        c = list(map(float, c))
        genes = [10, 20]#, 40, 50, 60]
        inds = [8]#, 6, 10]
        methods = [1,2]#, 4]
        probs = [1]#, 0.05, 0.2]
        species = [2, 4, 6]#, 3, 5]
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
        c = [1, 4, 8, 12, 16, 20]
        genes = [10, 20, 30]
        inds = [8]
        methods = [1, 2]
        probs = [4, 8, 16]
        species = [2, 4, 6, 8]
        trees = [3] # number of species trees
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
    batch_run = partial(setup, batch_folder, methods, probs, flow_dict, force)
    batch_iterator = product(depth, genes, inds, pop_size, species, trees)
    
    # choose run method
    run_parallel= lambda: tools.parmap(batch_run, batch_iterator)
    run_serial = lambda: list(map(batch_run, batch_iterator))
    f = run_parallel if parallel else run_serial

    # run experiment
    tools.timeit(f, "solving all problems")

