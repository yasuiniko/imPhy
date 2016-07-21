"""
Usage: batch.py all <batch_folder> <c> [-f] -g=n_gene_trees -i=n_ind -m=method... [-n=Ne] -p=prob... -s=n_sp -t=n_sp_trees
       batch.py generate <batch_folder> <c> [-f -g=n_gene_trees -i=n_ind -m=method... -n=Ne -p=prob... -s=n_sp --plus]
       batch.py drop <batch_folder> <n_sp> [-f -m=method... -p=prob --plus]
       batch.py impute <batch_folder> [-f -m=method... --plus]
       batch.py analyze <batch_folder> [-f]

Options:
  -h, --help            Show this help message.
  -f                    Perform operations in a way that overwrites
                        existing files.
  <batch_folder>        Path to exp_folder/batch_folder. Will contain
                        all information pertaining to the batch.
  --plus                Perform all steps following the given step. 
                        Example: batch.py drop my_batch 5 --plus
                            1) Drop leaves from trees in
                               my_batch/nexus/ and put their distance
                               matrices into my_batch/data/
                            2) Impute files in my_batch/data/
                            3) Analyze files in my_batch/solutions to 
                               create my_batch/stats

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
import itertools
import os
import shutil
import subprocess
import time

from analyze import analyze
from generateTrees import generateTrees
from tools import timeit

def generate(nexus, d, n_gene_trees, n_sp, n_ind, Ne, n_sp_trees):
    
    # call generate trees
    f = lambda: generateTrees([d], 
                              n_gene_trees,
                              n_sp,
                              n_ind,
                              Ne,
                              n_sp_trees,
                              nexus)
    timeit(f, "generating trees")

def drop(batch_folder, nexus, data, n_sp, prob_missing):
    
    # function to call RandomGenerator.R
    def call_R(args):
        p_drop, fname = args
        in_name = os.path.join(nexus, fname)
        namelist = fname[:-4].split("_")
        namelist.append("p{}".format(p_drop))
        out_name = os.path.join(data, "_".join(namelist))

        subprocess.check_call(["Rscript",
                               "RandomGenerator.R",
                               in_name,
                               "-o", out_name,
                               "-p{}".format(p_drop),
                               "-s{}".format(n_sp)])

    # identify names of gene trees
    gene_trees = lambda fname: '_e' in fname.lower()
    basenames = filter(gene_trees, os.listdir(nexus))

    # call RandomGenerator.R
    f = lambda: list(map(call_R,
                         itertools.product(prob_missing, basenames)))
    timeit(f, "dropping leaves")

def impute(batch_folder, data, solutions, batch_base, methods):

    # function to move or copy files
    def move(dest, s="move"):
        if s == "move":
            def to(src):
                os.rename(src, os.path.join(dest, os.path.basename(src)))
        else:
            def to(src):
                shutil.copy(src, dest)
        return to

    # function to call optimizer
    def call_imp(args):
        basename, method = args
        nameroot = basename[:-4]
        methodstr = "" if method == 1 else 2
        program = "./missing{}.o".format(methodstr)
        
        f = lambda: subprocess.check_call([program, nameroot])
        timeit(f, "imputing {}".format(nameroot))

        namelist = nameroot.split("_")
        namelist.insert(-1, "m{}".format(method))
        dest = "_".join(namelist) + ".sol"
        src = nameroot + ".sol"
        cpp_sol = os.path.abspath('sol')
        os.rename(os.path.join(cpp_sol, src), os.path.join(cpp_sol, dest))

    # set up
    try:
        os.chdir('cpp')
    except FileNotFoundError as e:
        print("Please ensure that you call this file from the src " +
               "folder, which contains the cpp folder.")
        raise e

    # set up files for optimizer
    cpp_data = os.path.abspath('data')
    cpp_sol = os.path.abspath('sol')
    dropped_dists = lambda fname: fname[-8:-4] != 'true'
    basenames = list(filter(dropped_dists, os.listdir(data)))
    in_basenames = lambda f: f[:-4] in map(lambda x: x[:-4], basenames)

    # copy files into cpp_data
    data_path = lambda f: os.path.join(data, f)
    list(map(move(cpp_data, "copy"), map(data_path, os.listdir(data))))

    # impute files
    f = lambda: list(map(call_imp, itertools.product(basenames, methods)))
    timeit(f, "imputing {} problems".format(len(basenames)))

    # delete files in cpp_data
    c_data_path = lambda f: os.path.join(cpp_data, f)
    in_data = lambda f: f[:-4] in map(lambda x: x[:-4], os.listdir(data))
    list(map(os.remove,
             map(c_data_path, 
                 filter(in_data,
                        os.listdir(cpp_data)))))

    # move files to solutions
    c_sol_path = lambda f: os.path.join(cpp_sol, f)
    in_sols = lambda f: set(f.split("_")[:4]) <= set(batch_base.split("_"))
    list(map(move(solutions, "move"),
             map(c_sol_path,
                 filter(in_sols,
                        os.listdir(cpp_sol)))))
    os.chdir('..')

def make_flow(flow_dict):
    flow = ["generate", "drop", "impute", "analyze"]
    if not flow_dict["all"]:
        while not (flow_dict[flow[0]] if flow else True): flow.pop(0)
        if not flow_dict['--plus']:
            flow = [flow[0]]

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
            drop(batch_folder, nexus, data, n_sp, prob_missing)

    # impute
    if current_step_is("impute"):
        if force or not exists(batch_folder, "solutions", m_tags):
            impute(batch_folder, data, solutions, batch_base, methods)

    # get summary stats for each file
    if current_step_is("analyze"):
        if force or not exists(batch_folder, "stats"):
            analyze(batch_folder)

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

    run_batch(batch_folder,
              species_depth,
              n_gene_trees,
              n_ind,
              n_sp,
              n_sp_trees,
              Ne,
              prob_missing,
              methods,
              args,
              force)
