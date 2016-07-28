"""
Usage: analyze.py <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders.
"""

import decimal
import dendropy # pip
import docopt # conda
from functools import partial, reduce
import io
import itertools
import numpy as np # conda
import math
import operator
import os
import pandas as pd # conda
import subprocess
import tempfile

import tools

# DendroPy naming
pdm = dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix

def rf_dist(args):
    assert len(args) == 2
    return dendropy.calculate.treecompare.symmetric_difference(*args)

def bhv_dist(args, batch_folder):
    assert len(args) == 2

    t2n = lambda tree: tree.as_string(schema='newick',
                                      suppress_leaf_taxon_labels=False,
                                      suppress_leaf_node_labels=True,
                                      suppress_internal_taxon_labels=True,
                                      suppress_internal_node_labels=True,
                                      suppress_rooting=True)

    (fd, infile) = tempfile.mkstemp()
    outfile = infile + 'out'

    try:
        tfile = os.fdopen(fd, "w")
        list(map(tfile.write, map(t2n, args)))
        tfile.close()
        
        gtp = "java^^-jar^^gtp.jar^^-o^^{}^^{}".format(outfile, infile)
        subprocess.check_call(gtp.split("^^"))

    finally:
        os.remove(infile)

    with open(outfile, 'r') as f:
        distance = to_decimal(f.readline().split()[-1])

    os.remove(outfile)

    return distance

def to_decimal(x):
    return decimal.Decimal(str(x))

def vector2list(v):
    """
    Converts the vectorized form of a distance matrix (in the form of 
    a numpy array) to a flattened distance matrix (in the form of a
    list) filled with Decimals. 
    """

    # get dimensions of matrix
    n = int((1 + math.sqrt(8 * len(v) + 1)) / 2)
    m = np.zeros((n, n))

    # fill upper triangle
    m[np.triu_indices(n, k=1)] = v

    # fill lower triangle, convert to list, convert to elements to Decimal
    return list(map(lambda x: decimal.Decimal(str(x)), (m + m.T).flatten()))

def make_taxa(nexus_file):
    taxa = []
    with open(nexus_file, 'r') as f:
        read = False
        for line in f:
            if "TAXLABELS" in line:
                read = True
            elif ";" in line:
                read = False
            elif read:
                taxa.append(line.strip())

    return dendropy.TaxonNamespace(taxa, label=nexus_file)

def reconstruct_trees(taxa, vectors):
    """
    Creates trees using DendroPy's Neighbor Joining and UPGMA methods.
    """

    def helper(vector):
        csv_string = ",{}\n".format(",".join(taxa.labels()))

        # get dimensions of matrix
        n = int((1 + math.sqrt(8 * len(vector) + 1)) / 2)
        m = np.zeros((n, n))

        # fill upper triangle
        m[np.triu_indices(n, k=1)] = vector

        # fill csv_string
        label = iter(taxa.labels())
        for row in m:
            csv_string += next(label) + ','
            for i, cell in enumerate(row):
                csv_string += str(cell)
                if i < len(row) - 1:
                    csv_string += ','
            csv_string += "\n"
        
        # get distance matrix
        dist_matrix = pdm.from_csv(io.StringIO(csv_string), 
                                   taxon_namespace=taxa,
                                   delimiter=',',
                                   is_allow_new_taxa=False)
        return dist_matrix

    # make trees
    nj = lambda x: x.nj_tree()
    upgma = lambda x: x.upgma_tree()
    distances = list(map(helper, vectors))

    return map(nj, distances), map(upgma, distances)

def percentiles_of(x):
    """
    Returns a precentile function for list x. 
    """
    x = sorted(x)
    n = len(x)
    
    def percentile(p):
        """
        Calculates the pth percentile of x.
        """
        p = p/100
        if 0<=p<=1/(n+1):
            index = 0
        elif 1/(n+1) < p < n/(n+1):
            index = p*(n+1) - 1
        else:
            index = n - 1

        ind = math.floor(index)
        index = decimal.Decimal(str(index))

        return x[ind] + index%1 * (x[ind+1] - x[ind]) if index % 1 else x[ind]

    return percentile

def calc(error, mode="stats", sol_file=""):
    # calculations
    sq_err = list(map(lambda x: to_decimal(x)**2, error))
    imp_sq_err = list(filter(lambda x: x != 0, sq_err))
    imp_err = list(map(lambda x: x.sqrt(), imp_sq_err))
    sse = to_decimal(sum(imp_sq_err))
    n = len(imp_sq_err) if mode == "stats" else len(sq_err)
    mse = sse/n if n else 0
    rmse = math.sqrt(mse)

    retval = None

    if mode == "stats":
        # calculate normalized error
        def param_value(tag):
            return list(filter(lambda tags: tag in tags,
                               os.path.basename(sol_file).split("_")))[-1][1:]
        species_depth = int(param_value('d'))
        theoretical_max = 2 * species_depth 
        nrmse = rmse / theoretical_max

        # calculate percentiles
        p = percentiles_of(imp_err)

        retval = [p(0), p(25), p(50), p(75), p(100), rmse, nrmse]

    elif mode == 'rf': 
        retval = rmse / (2 * (n - 3))

    elif mode == 'bhv':
        retval = rmse

    return retval

def get_stats(solution_file, true_file):
    """
    Finds quartiles, Root Mean Squared Error (RMSE), and Normalized
    Root Mean Squared Error (NRMSE) using files with solution vectors
    and actual distance vectors.
    """
    # read file
    sol_list = np.loadtxt(solution_file)
    true_list = np.loadtxt(true_file)

    if sol_list.size == 0:
        return ["Solution file was empty. Please check for imputation errors."]

    # get lists of distances
    imp_dists = list(itertools.chain.from_iterable(map(vector2list, sol_list)))
    og_dists = list(itertools.chain.from_iterable(map(vector2list, true_list)))

    # calculations
    err = map(lambda x: x[0] - x[1], zip(imp_dists, og_dists))
    return calc(err, sol_file=solution_file)

def tree_stats(solution_file, true_file, tree_gen, batch_folder):    
    # read file
    sol_list = np.loadtxt(solution_file)
    true_list = np.loadtxt(true_file)

    # get trees
    imp_nj, imp_upgma = tree_gen(sol_list)
    og_nj, og_upgma = tree_gen(true_list)
    nj_trees = list(zip(imp_nj, og_nj))
    upgma_trees = list(zip(imp_upgma, og_upgma))

    # get distances
    bhv = partial(bhv_dist, batch_folder=batch_folder)
    bhv_nj = calc(map(bhv, nj_trees), mode='bhv')
    rf_nj = calc(map(rf_dist, nj_trees), mode='rf')
    bhv_upgma = calc(map(bhv, upgma_trees), mode='bhv')
    rf_upgma = calc(map(rf_dist, upgma_trees), mode='rf')

    return bhv_nj, rf_nj, bhv_upgma, rf_upgma
    # return rf_nj_t, rf_upgma_t

def write_to(batch_folder, solution_file, desc=""):
    """
    Creates a function that writes x to the relevant file in the stats
    subdirectory. 
    """
    stats = os.path.join(batch_folder, "stats")
    if not os.path.isdir(stats):
        os.makedirs(stats)

    basename = os.path.basename(solution_file)[:-4] + desc

    def write(x):
        """
        Writes x to file.
        """
        if type(x) == tuple or type(x) == list:
            x = list(itertools.chain(x))
        with open(os.path.join(stats, basename + '.txt'), 'w') as f:
            f.write(" ".join(map(str, x)))
    return write

def match(infile, file_list, tag):
    """
    Finds the first file in file_list containing the same tag value as 
    infile.
    """
    get_tag = lambda f: list(filter(lambda x: tag in x, f[:-4].split("_")))[-1]
    return list(filter(lambda f: get_tag(infile) in f, file_list))[0]

def analyze(batch_folder):
    """
    File to write a file in the stats subfolder for each file in the
    solutions file.
    """
    # functions to get files
    def get_files(subdir):
        abs_path = lambda fname: os.path.join(batch_folder, subdir, fname)
        subdir_path = os.path.join(batch_folder, subdir)
        return map(abs_path, os.listdir(subdir_path))
    
    def tagged_files(subdir, tag):
        has_tag = lambda fpath: tag in os.path.basename(fpath)
        return sorted(list(filter(has_tag, get_files(subdir))))

    # get files 
    true_files = tagged_files("data", "_true.txt")
    sol_files = tagged_files("solutions", ".sol")
    nexus_files = tagged_files("nexus", "_e")

    # match solution files to files from which they were created
    match_sol_true = lambda sol: match(sol, true_files, 'p')

    # match solution files with a tree generation function
    sol2nexus = {sol: match(sol, nexus_files, 'e') for sol in sol_files}
    nex2trees = {n: partial(reconstruct_trees, make_taxa(n)) for n in nexus_files}
    match_sol_treegen = lambda sol: nex2trees[sol2nexus[sol]]

    # make function to write stats
    def write_stats(tup): 
        write_to(batch_folder, tup[0])(get_stats(*tup))
    def write_trees(tup):
        write_to(batch_folder, tup[0], '_tree')(tree_stats(*tup, batch_folder))

    # write stats
    stats_files = zip(sol_files, map(match_sol_true, sol_files))
    tree_files = zip(sol_files,
                     map(match_sol_true, sol_files), 
                     map(match_sol_treegen, sol_files))
    list(map(write_stats, stats_files))
    list(map(write_trees, tree_files))

def summary(exp_folder):
    """
    Creates a csv file with summary statistics from each batch_folder in 
    exp_folder.
    """

    # functions to parse the exp_folder
    def values_from_tags(tags):
        """
        Returns a list of lists of the values of each parameter in tags.
        """

        # check if string is a valid parameter
        valid_param = lambda s: s[0].isalpha() and s[1].isdigit()
        # split by underscore
        def split(s):
            if list(filter(lambda x: x.isdigit(), s[-4:])):
                return list(filter(lambda x: x, s.split("_")))
            else:
                return list(filter(lambda x: x, s[:-4].split("_")))

        param_list = set(filter(valid_param,
                                itertools.chain.from_iterable(map(split,
                                                                  tags))))  

        def values_from_tag(param_tag):
            """
            Get all values of a parameter in param_list that is
            uniquely identified by its param_tag. 
            """
            # get the value of a parameter
            param_value = lambda param: param.split(param_tag)[-1]
            # check if the parameter has the correct tag
            correct_tag = lambda param: param_tag in param
            
            return list(map(param_value, 
                            filter(correct_tag, 
                                   param_list)))
        
        # get the first element of a string
        first_char = lambda s: s[0]

        # sorted list of parameter tags
        param_tags = sorted(list(set(map(first_char, param_list))))

        # list of lists of values of parameters
        return list(map(values_from_tag, param_tags))

    def get_row(vect, param_values, rowsize, modifier=""):
        """
        Get the row of values to store in the CSV file according to
        vect, a vector describing the values of parameters in the
        experiment.
        """

        # get values of current parameters
        get_vals = lambda x: param_values[x[0]][x[1]]
        params = list(map(get_vals, enumerate(vect)))

        # get paths
        batch_folder = tools.batch.format(*params)
        nameroot = tools.fileroot.format(*params)
        filename = nameroot + modifier + ".txt"
        stats_path = os.path.join(exp_folder,
                      batch_folder,
                      "stats",
                      filename)

        # get data
        try:
            # stats_path points to non-tree stats file
            row = params + list(np.loadtxt(stats_path))

        except ValueError as e:
            # stats_path points to tree stats file or imputation
            # error
            with open(stats_path, 'r') as f:
                cln = lambda x: x.translate({ord(c): None for c in '[],'})
                data = list(map(cln, f.readline().split()))
            
            if list(filter(lambda x: x.isdigit(), data)):
                row = params + data
            else:
                row = [float('nan')] * rowsize

        return row

    # make sure exp_folder is absolute path
    exp_folder = os.path.abspath(exp_folder)

    # specify column names (for later use)
    cols = ["Species Depth", "Species Tree Index", "Number of Gene Trees", 
            "Number of Individuals per Species", "Method",
            "Effective Population Size", "Leaf Dropping Probability",
            "Number of Species", "Abs Imputation Error (AIE) min", 
            "AIE lower quartile", "AIE median", "AIE upper quartile",
            "AIE max", "Imputation RMSE", "Imputation NRMSE"]

    # figure out which experiments were run
    valid = lambda name: name[0].isalpha() and name[1].isdigit()
    supertags = list(filter(valid, os.listdir(exp_folder)))
    subtags = list(filter(valid, os.listdir(os.path.join(exp_folder, 
                                                         supertags[0],
                                                         'stats'))))

    split_tags = itertools.chain.from_iterable(map(lambda x: x.split("_"), 
                                                   supertags))
    n_gene_trees = max(list(map(lambda x: int(x[1:]),
                                filter(lambda x: 'g' in x, 
                                       split_tags))))

    tree_cols = ["Species Depth", "Species Tree Index", "Number of Gene Trees", 
            "Number of Individuals per Species", "Method",
            "Effective Population Size", "Leaf Dropping Probability",
            "Number of Species", "bhv_nj", "rf_nj", "bhv_upgma", "rf_upgma"]
    # tree_types = ["bhv_nj", "rf_nj", "bhv_upgma", "rf_upgma"]
    # lst_sum = lambda a, b: a + b
    # get_names = lambda x: [x]*n_gene_trees
    # tree_cols = base_cols + reduce(lst_sum, map(get_names, tree_types), [])

    # get values of parameters used in this experiment
    param_values = values_from_tags(supertags + subtags)
    coord_data = partial(get_row,
                         param_values=param_values,
                         rowsize=len(cols),
                         modifier="")
    coord_tree = partial(get_row, 
                         param_values=param_values,
                         rowsize=len(tree_cols), 
                         modifier="_tree")

    # get sizes to create the numpy array
    dimensions = list(map(len, param_values))
    size = reduce(operator.mul, dimensions, 1)

    # build data arrays
    data = np.empty((size, len(cols)), dtype=float)
    tree = np.empty((size, len(tree_cols)), dtype=float)
    coordinates = itertools.product(*list(map(range, dimensions)))

    for ind, coordinate in enumerate(coordinates):
        data[ind] = coord_data(coordinate)
        tree[ind] = coord_tree(coordinate)

    # write data
    outpath = os.path.join(exp_folder, "summary.csv")
    pd.DataFrame(data, columns=cols).to_csv(outpath)
    # write trees
    treepath = os.path.join(exp_folder, "tree_summary.csv")
    pd.DataFrame(tree, columns=tree_cols).to_csv(treepath)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    infolder = args['<infolder>']

    summary(infolder)