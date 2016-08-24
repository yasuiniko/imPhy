"""
Usage: experiment.py <exp_folder> [options]

Options:
  <exp_folder>          Path to destination folder for the sets of data.

  -c                    Compiles statistics about the experiment after
                        running. Requires matplotlib, pandas, and
                        seaborn.  
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
    comp_stats = args['-c']
    force = args['-f']
    parallel = args['-p']
    test, experiment = args['-t'], not args['-t']

    # set up filesystem
    if exp_folder[0] != '/':
        exp_folder = os.path.abspath(os.path.join("..", exp_folder))
    if not os.path.isdir(exp_folder):
        os.makedirs(exp_folder)
    heatpath = os.path.join(exp_folder, 'heatmaps')
    if os.path.isdir(heatpath):
        shutil.rmtree(heatpath)
    os.makedirs(heatpath)

    # set up logger
    logpath = os.path.join(exp_folder, "output.log")
    logger = tools.MultiLogger(logpath)

    # test run
    if test:
        c = [20]#, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
        c = list(map(float, c))
        genes = [3,4,5]#, 40, 50, 60]
        inds = [8]#, 6, 10]
        methods = [1]#, 4]
        probs = [8]#, 0.05, 0.2]
        species = [2]#, 3, 5]
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
        c = [1, 4, 8, 12, 16, 20]  # c ratio
        genes = [200, 500]         # number of genes
        inds = [8]                 # number of individuals per species
        methods = [1, 2]           # imputation methods to use
        probs = [4, 8, 16]         # leaf dropping probabilities/denominators
        species = [2, 4, 6, 8]     # number of species 
        trees = [3]                # number of species trees
        pop_size = [10000]         # effective population size

        # convert c list from integers to floats
        c = list(map(float, c))
        # depth (in generations) is equal to c*pop_size, for each 
        # combination of values of c and pop_size. If you prefer, 
        # replace the following line with a hard-coded depth list, as 
        # seen above in the definition for c. 
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
    def run_parallel():
        tools.parmap(batch_run, batch_iterator)

    def run_serial():
        for batch in batch_iterator:
            batch_run(batch)

    f = run_parallel if parallel else run_serial

    # run experiment
    tools.timeit(f, "solving all problems", logger.getLogger(__name__))

    if comp_stats:
        compile_stats(exp_folder)
