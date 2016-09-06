"""
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: experiment.py <exp_folder> [options]

Options:
  <exp_folder>          Path to destination folder for the sets of data.

  -d                    Creates diagnostic plots about the experiment 
                        after running.                  [default: False]
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

from compile_stats import compile_stats, make_plots
from batch import run_batch
import settings
import tools

def setup(batch_folder, methods, probs, flow_dict, dists, force, tup):

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
              dists=dists,
              flow_dict=flow_dict,
              force=force)

if __name__ == "__main__":
    args = docopt.docopt(__doc__)

    # get args
    exp_folder = args['<exp_folder>']
    diag_plots = args['-d']
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
        opts = settings.Settings.get_t()
    # experimental set up
    if experiment:
        opts = settings.Settings.get_e()
    
    # unpack opts
    methods = opts.methods
    probs = opts.probs
    flow_dict = opts.flow_dict
    dists = opts.dists
    depths = opts.depths
    genes = opts.genes
    inds = opts.inds
    pop_sizes = opts.pop_sizes
    species = opts.species
    trees = opts.trees

    # check java installation
    needs_java = flow_dict["analyze"]
    if needs_java:
        tools.get_output(['java', '-version'])

    # batch folder naming scheme
    batch_folder = os.path.join(exp_folder, tools.batch_general)
    setup_args = [batch_folder, methods, probs, flow_dict, dists, force]
    batch_run = partial(setup, *setup_args)
    batch_iterator = product(depths, genes, inds, pop_sizes, species, trees)
    
    # choose run method
    def run_parallel():
        tools.parmap(batch_run, batch_iterator)

    def run_serial():
        for batch in batch_iterator:
            batch_run(batch)

    f = run_parallel if parallel else run_serial

    # run experiment
    tools.timeit(f, "solving all problems", logger.getLogger(__name__))
    compile_stats(exp_folder, dists)

    if diag_plots:
        make_plots(exp_folder)
