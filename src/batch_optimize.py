"""
Usage: batch_optimize.py generate <outfolder> <c>... [-g=n_gene_trees -i=n_ind -m=method -n=Ne -p=prob -s=n_sp --simple]
       batch_optimize.py drop <infolder> <n_sp> [-m=method -p=prob --simple]
       batch_optimize.py impute [-m=method -t=tag]
       batch_optimize.py impute_all [-m=method]

Options:
  -h, --help            Show this help message.
  --simple              Only carry out the given command (no imputing).

  Tree Generation:
  <outfolder>           Destination folder for generated trees. 
  <c>                   Values of SD:Ne to consider.
  -g, --n_gene_trees=.  Number of gene trees.            [default: 1000]
  -i, --n_ind=.         Number of individuals per species.  [default: 2]
  -n, --Ne=.            Effective population size.      [default: 10000]
  -s, --n_sp=.          Number of species.                  [default: 5]

  Leaf Dropping:
  <infolder>            Folder from which to source the gene trees.
  <n_sp>                Number of species.
  -p, --prob=.          Probability of missingness/leaf.  [default: 0.2]

  Imputation:
  -m, --method=.        Method to impute files (1 or 2).    [default: 1]
  -t, --tag=.           All files with this tag will be imputed.
"""

import os
import shutil
import subprocess
import docopt
import time
import itertools

from generateTrees import generateTrees

def timeit(f, s):
    bigs = s[0].upper() + s[1:]
    smalls = s[0].lower() + s[1:]
    print("{}...".format(bigs))
    t = time.time()
    f()
    print("Done {} in {}s.\n".format(smalls, time.time() - t))

def get_data(fname):
    return fname[-8:-4] != 'true'

def get_gene_trees(fname):
    valid = True
    if 'gene' not in fname and "Gene" not in fname:
        # print("Skipping file {} since it".format(fname) + 
        #       " does not contain 'gene'.")
        valid = False
    return valid

def drop_leaves(n_sp, infolder, prob_missing):

    def helper(fname):
        in_name = os.path.join(infolder, fname)
        out_name = os.path.join(os.path.abspath('cpp/data'), 
                                "{}_{}".format(os.path.basename(infolder),
                                               fname[:-4]))
        subprocess.check_call(["Rscript",
                               "RandomGenerator.R",
                               in_name,
                               "-o", out_name,
                               "-p{}".format(prob_missing),
                               "-s{}".format(n_sp)])
    return helper

def impute(method, prefix):
    assert method == 1 or method == 2
    method = "" if method == 1 else method

    def helper(fname):
        name = "{}{}".format(prefix, fname[:-4])
        program = "./missing{}.o".format(method)
        
        f = lambda: subprocess.check_call([program, name])
        timeit(f, "imputing {}".format(name))
    return helper

def cleanup(tag_names, infolder, prefix):

    abs_join = lambda *args: os.path.abspath(os.path.join(*args))

    # make subdirectories
    nex = abs_join(infolder, "nexus")
    os.makedirs(nex)
    data = abs_join(infolder, "data")
    os.makedirs(data)
    sol = abs_join(infolder, "solutions")
    os.makedirs(sol)

    # move nex files
    list(map(lambda name: shutil.copy(os.path.abspath(name), nex),
             filter(os.path.isfile, os.listdir(infolder))))

    # move data files
    missing = lambda fname: abs_join("data", 
                                     "{}{}.txt".format(prefix, fname[:-4]))
    actual = lambda fname: abs_join("data", 
                                    "{}{}_true.txt".format(prefix, fname[:-4]))
    data_files = map(lambda fname: [missing(fname), actual(fname)], tag_names)
    list(map(lambda name: shutil.copy(name, data),
             itertools.chain.from_iterable(data_files)))

    # move solution files
    list(map(lambda name: shutil.copy(name, sol),
             map(lambda fname: abs_join("sol", "{}{}.sol".format(prefix, 
                                                                 fname[:-4])),
                 filter(lambda x: 'true' not in x, tag_names))))

if __name__ == "__main__":
    # get args
    args = docopt.docopt(__doc__)

    try:
        outfolder = os.path.abspath(args['<outfolder>'])
    except:
        outfolder = None
    try:
        infolder = os.path.abspath(args['<infolder>'])
    except:
        infolder = None
    c = list(map(float, args['<c>']))
    n_gene_trees = int(args['--n_gene_trees'])
    n_ind = int(args['--n_ind'])
    n_sp = int(args['--n_sp']) if args['<n_sp>'] is None else args['<n_sp>']
    prob_missing = float(args['--prob'])
    Ne = int(args['--Ne'])
    method = int(args['--method'])
    tag = args['--tag']

    # generate trees
    if args['generate']:
        
        f = lambda: generateTrees(c, n_gene_trees, n_sp, n_ind, Ne, outfolder)
        timeit(f, "generating trees")

    # drop leaves
    if args['drop'] or (args['generate'] and not args['--simple']):

        # setup
        infolder = infolder if outfolder is None else outfolder
        basenames = filter(get_gene_trees, os.listdir(infolder))

        # drop leaves
        f = lambda: list(map(drop_leaves(n_sp, infolder, prob_missing), 
                             basenames))
        timeit(f, "dropping leaves")


    # impute by tag or data from earlier in the session
    if not args['impute_all'] and not args['--simple']:

        # set up
        os.chdir('cpp')
        tag = '' if tag is None else tag
        if infolder is None:
            infolder = os.path.abspath('cpp/data')
            basenames = filter(get_data, os.listdir('data'))
        else:
            basenames = filter(get_gene_trees,
                               os.listdir(infolder))
        tag_names = list(filter(lambda x: tag in x, basenames))
        prefix = os.path.basename(infolder)
        prefix = "" if prefix == 'data' else "{}_".format(prefix)

        # impute based on tag
        f = lambda: list(map(impute(method, prefix), tag_names))
        timeit(f, "imputing {} problems".format(len(tag_names)))

        if not args['impute']:
            f = lambda: cleanup(tag_names, infolder, prefix)
            # timeit(f, "cleaning up")

    # imput all
    if args['impute_all']:

        # set up
        os.chdir('cpp')
        
        # impute all

        f = lambda: list(map(impute(method, 'all'),
                             filter(get_data, os.listdir('data'))))
        timeit(f, "imputing all")
