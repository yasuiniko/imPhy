"""
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: batch.py all <batch_folder> -d=dists... <c> [-f] -g=n_gene_trees -i=n_ind -m=method... [-n=Ne] -p=prob... -s=n_sp -t=n_sp_trees
       batch.py generate <batch_folder> -d=dists... <c> [-f -g=n_gene_trees -i=n_ind -m=method... -n=Ne -p=prob... -s=n_sp]
       batch.py drop <batch_folder> -d=dists... <n_sp> [-f -m=method... -p=prob]
       batch.py impute <batch_folder> -d=dists... [-f -m=method...]
       batch.py analyze <batch_folder> -d=dists... [-f]

Options:
  -h, --help            Show this help message.
  -f                    Perform operations in a way that overwrites
                        existing files.
  <batch_folder>        Path to exp_folder/batch_folder. Will contain
                        all information pertaining to the batch.
  -d, --dists=.         Distance metrics to use (-d bhv -d norm -d rf)

  Tree Generation:
  <d>                   Values of Species depth to consider.
  -g, --n_gene_trees=.  Number of gene trees.            [default: 1000]
  -i, --n_ind=.         Number of individuals per species.  [default: 2]
  -n, --Ne=.            Effective population size.      [default: 10000]
  -s, --n_sp=.          Number of species.                  [default: 5]
  -t, --n_sp_trees=.    Number of species trees.            [default: 2]

  Leaf Dropping:
  <n_sp>                Number of species.
  -p, --prob=.          Probability of missingness/leaf.  [default: 0.2]

  Imputation:
  -m, --method=.        Method to impute files (1 or 2).    [default: 1]
"""

import docopt
from functools import partial
import gc
import itertools
import logging
import os
import shutil
from subprocess import Popen, PIPE, CalledProcessError
import time

from compile_stats import analyze
from generateTrees import generateTrees
from tools import ext, ext_len, gzip_to, gunzip_to, timeit, get_output

def generate(nexus, d, n_gene_trees, n_sp, n_ind, Ne, n_sp_trees):
    
    # call generate trees
    f = lambda: generateTrees([d], 
                              n_gene_trees,
                              n_sp,
                              n_ind,
                              Ne,
                              n_sp_trees,
                              nexus)
    timeit(f, "generating trees", logging.getLogger("generate_trees"))

def drop(batch_folder, nexus, data, n_sp, prob_missing, force):
    
    # function to call RandomGenerator.R
    def call_R(args):

        # setup
        p_drop, fname = args
        in_name = os.path.join(nexus, fname)
        namelist = fname[:-4].split("_")
        namelist.append("p{}".format(p_drop))
        out_name = os.path.join(data, "_".join(namelist))

        # remove old files if force is true
        if force:
            if os.path.exists(out_name+".txt"):
                os.remove(out_name+".txt")
            if os.path.exists(out_name+"_true.txt"):
                os.remove(out_name+"_true.txt")

        # call the R script               
        get_output(["Rscript",
                    "RandomGenerator.R",
                    in_name,
                    "-o", out_name,
                    "-p{}".format(p_drop),
                    "-s{}".format(n_sp)],
                    logging.getLogger("drop"))

        # compress data files
        dest_root = os.path.split(out_name)[0]
        gzip_to(dest_root, out_name+'.txt')
        gzip_to(dest_root, out_name+'_true.txt')

    # identify names of gene trees
    gene_trees = lambda fname: '_e' in fname.lower()
    basenames = filter(gene_trees, os.listdir(nexus))

    # call RandomGenerator.R
    f = lambda: list(map(call_R, itertools.product(prob_missing, basenames)))
    timeit(f, "dropping leaves", logging.getLogger("drop"))

def impute(batch_folder, data, solutions, batch_base, methods):
    """
    Imputes the batch.
    """

    def call_imp(args):
        """
        Calls the cpp imputation software.
        """
        
        # set up variables
        basename, method = args
        nameroot = basename[:-ext_len]
        program = "./missing{}.o".format(method)
        
        # call the imputation software
        
        def imp_wrapper():
            """
            Wraps the imputation software call in a try/catch to
            continue imputing other files in the batch if some calls 
            fail.
            """
            logger = logging.getLogger("impute")
            try:
                p = Popen([program, nameroot], 
                                     stdout=PIPE,
                                     stderr=PIPE)
                output, err = p.communicate()
                
                if output:
                    logger.debug("Output of {}: {}".format([program, nameroot],
                                                           output))

            except CalledProcessError:
                p.terminate()
                logger.error("Imputation error for basename '{}'".format(args[0]) + 
                         " using method '{}'.".format(args[1]))

            except OSError as e:
                logger.error("OSError while imputing {}".format(nameroot))
                raise e

        timeit(imp_wrapper,
               "imputing {}".format(nameroot),
               logging.getLogger("impute"))

        # garbage collection
        gc.collect()

        # rename solution file to include the method number
        namelist = nameroot.split("_")
        namelist.insert(-1, "m{}".format(method))
        dest = "_".join(namelist) + ".sol"
        src = nameroot + ".sol"
        cpp_sol = os.path.abspath('sol')

        logger = logging.getLogger("impute")
        try:
            os.rename(os.path.join(cpp_sol, src), os.path.join(cpp_sol, dest))
            logger.debug("Successfully renamed solution file for basename '{}' and method '{}'".format(*args))
        except FileNotFoundError as e:
            raise e

    # set up
    try:
        os.chdir('cpp')
    except FileNotFoundError as e:
        logger = logging.getLogger("impute")
        logger.error("Please ensure that you call this file from the src " +
               "folder, which contains the cpp folder.")
        raise e

    # set up files for optimizer
    cpp_data = os.path.abspath('data')
    cpp_sol = os.path.abspath('sol')
    dropped_dists = lambda fname: '_true.' not in fname and "txt.gz" in fname
    basenames = list(filter(dropped_dists, os.listdir(data)))

    # copy files into cpp_data
    add_data_path = lambda f: os.path.join(data, f)
    data_paths = map(add_data_path, basenames)
    to_cpp_data = partial(gunzip_to, cpp_data)
    list(map(to_cpp_data, data_paths))

    # impute files
    logger = logging.getLogger("impute")
    for args in itertools.product(basenames, methods):
        logger.debug("Imputing {} with method {}".format(*args))
        call_imp(args)

    # delete files in cpp_data
    c_data_path = lambda f: os.path.join(cpp_data, f)
    data_paths = list(map(lambda x: x[:-ext_len], os.listdir(data)))
    in_data = lambda f: f[:-4] in data_paths
    list(map(os.remove,
             map(c_data_path, 
                 filter(in_data,
                        os.listdir(cpp_data)))))

    # move files to solutions
    c_sol_path = lambda f: os.path.join(cpp_sol, f)
    batch_ids = batch_base.split("_")
    in_batch = lambda f: set(f.split("_")[:len(batch_ids)]) <= set(batch_ids)
    correct_files = filter(in_batch, os.listdir(cpp_sol))
    correct_paths = map(c_sol_path, correct_files)
    to_sols = partial(gzip_to, solutions)
    list(map(to_sols, correct_paths))

    # move back up one level
    os.chdir('..')

def make_flow(flow_dict):
    """
    Finds the True entries in the flow_dict and returns a function that
    takes a step name and returns True if that step is True in flow_dict
    or None otherwise.
    """
    steps = ["generate", "drop", "impute", "analyze"]
    flow = list(filter(lambda step: flow_dict[step], steps))

    return lambda s: flow.pop(0) if flow and flow[0] == s else None

def exists(batch_folder, subfolder, tags=[""]):
    # batch folder and subfolder exist
    folders_exist = os.path.exists(os.path.join(batch_folder, subfolder))
    # subfolder is not empty
    non_empty = os.listdir(os.path.join(batch_folder, subfolder))
    
    # make sure all tags are present in the filenames of the subfolder
    for f in non_empty:
        for tag in tags:
            if tag in f:
                tags.remove(tag)
                break
    all_tags_found = True if tags == [] else False

    return True if folders_exist and non_empty and all_tags_found else False

def run_batch(batch_folder,
              species_depth,
              n_gene_trees,
              n_ind,
              n_sp,
              n_sp_trees,
              Ne,
              prob_missing,
              methods,
              flow_dict,
              dists,
              force=False):
    # set up folders
    directories = ["data", "nexus", "solutions", "stats"]
    dpath = lambda s: os.path.join(batch_folder, s)
    dpaths = list(map(dpath, directories))
    data, nexus, solutions, stats = dpaths
    batch_base = os.path.basename(batch_folder)

    # create folders if they don't already exist
    mdir = lambda d: None if os.path.isdir(d) else os.makedirs(d)
    list(map(mdir, dpaths))

    current_step_is = make_flow(flow_dict)

    tag = lambda s, lst: list(map(lambda x: s+str(x), lst))
    m_tags = tag("m", methods)
    p_tags = tag("p", prob_missing)

    # generate trees
    if current_step_is("generate"):
        if force or not exists(batch_folder, "nexus"):
            generate(nexus,
                     species_depth,
                     n_gene_trees,
                     n_sp,
                     n_ind,
                     Ne,
                     n_sp_trees)

    # drop leaves
    if current_step_is("drop"):
        if force or not exists(batch_folder, "data", p_tags):
            drop(batch_folder, nexus, data, n_sp, prob_missing, force)

    # impute
    if current_step_is("impute"):
        if force or not exists(batch_folder, "solutions", m_tags):
            impute(batch_folder, data, solutions, batch_base, methods)

    # get summary stats for each file
    if current_step_is("analyze"):

        if force or not exists(batch_folder, "stats", m_tags + p_tags):
            timeit(lambda: analyze(batch_folder, dists),
                   "analyzing {}".format(os.path.basename(batch_folder)),
                   logging.getLogger("analyze"))

if __name__ == "__main__":
    # get args
    args = docopt.docopt(__doc__)

    c = float(args['<c>']) if args['<c>'] else None
    n_gene_trees = int(args['--n_gene_trees'])
    n_ind = int(args['--n_ind'])
    n_sp = int(args['<n_sp>']) if args['<n_sp>'] else int(args['--n_sp'])
    n_sp_trees = int(args['--n_sp_trees'])
    prob_missing = list(map(float, args['--prob']))
    Ne = int(args['--Ne'])
    methods = list(map(int, args['--method']))
    batch_folder = os.path.abspath(args['<batch_folder>'])
    force = args['-f']
    species_depth = c * Ne

    dists = args['--dists']

    flow_dict = {"generate": True if args['all'] or args['generate'] else False,
                 "drop":True if args['all'] or args['drop'] else False,
                 "impute":True if args['all'] or args['impute'] else False,
                 "analyze":True if args['all'] or args['analyze'] else False}
    
    run_batch(batch_folder,
              species_depth,
              n_gene_trees,
              n_ind,
              n_sp,
              n_sp_trees,
              Ne,
              prob_missing,
              methods,
              flow_dict,
              dists,
              force)
