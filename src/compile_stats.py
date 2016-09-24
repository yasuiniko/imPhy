"""

This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: compile_stats.py <exp_folder> <dists>...

Options:
  -h, --help            Show this help message.
  <exp_folder>          Folder containing data, nexus, and solutions
                        subfolders.
  <dists>               Distances calculated.   [default: bhv, rf, norm]

Description:
Contains code for the 'analyze' step, as well as to put all the
statistics created during 'analyze' into "interleaf_error.csv" and
"intertree_error.csv" for easy graphing.
"""

import dendropy # pip
import docopt # conda
from functools import partial, reduce
import io
import itertools
import numpy as np # conda
import math
import matplotlib.pyplot as plt # conda
import operator
import os
import pandas as pd # conda
import random
import seaborn as sns # conda
import subprocess
import tempfile

from tools import *

# DendroPy naming
pdm = dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix

def l2_norm(trees):
    """
    Calculates the normalized l2 norm between two trees.
    """
    vectorize = lambda tree: [edge.length for edge in tree.edges()]

    clean_none = lambda length: length if length else 0
    clean = lambda vect: np.array(list(map(clean_none, vect)))

    norm = lambda v1, v2: np.linalg.norm(v1 - v2)

    return norm(*list(map(clean, map(vectorize, trees))))

def rf_dist(trees):
    """
    Calculates the normalized Robinson-Foulds distance between two trees.
    """
    assert len(trees) == 2
    dist = dendropy.calculate.treecompare.symmetric_difference(*trees)
    n = len(trees[0].leaf_nodes())

    return dist/(2*(n-3))

def bhv_dist(trees, batch_folder, rooted):
    """
    Wrapper for Owen and Provan's GTP code to calculate the normalized 
    BHV distance between two trees.
    """
    # make temporary outfiles for GTP
    (fd, infile) = tempfile.mkstemp()
    outfile = infile + 'out'

    try:
        # open temporary file to hold trees
        tfile = os.fdopen(fd, "w")
        list(map(tfile.write, map(newick, trees)))
        tfile.close()
        
        # call GTP
        root_arg = "" if rooted else "-u^^"
        gtp = "java^^-jar^^gtp.jar^^-v^^-n^^{}-o^^{}^^{}".format(root_arg,
                                                                 outfile, 
                                                                 infile)
        null_logger = logging.getLogger("null")
        null_logger.propagate = False
        output = get_output(gtp.split("^^"), null_logger)

    finally:
        os.remove(infile)

    # calculate distances
    codim, distance = -1, -1.0
    with io.StringIO(output) as f:
        for line in f:
            if "Combinatorial type" in line:
                codim = codimension(line.split()[-1])
            if "Geodesic distance between start and target tree is" in line:
                distance = to_decimal(line.split()[-1])

    # remove outfile
    os.remove(outfile)

    return [distance, codim]

def newick(tree):
    return tree.as_string(schema='newick',
                          suppress_leaf_taxon_labels=False,
                          suppress_leaf_node_labels=True,
                          suppress_internal_taxon_labels=True,
                          suppress_internal_node_labels=True,
                          suppress_rooting=True)

def codimension(comb_block):
    """
    Find longest subblock in 'combinatorial types' block
    """
    # remove unnecessary characters
    clean = str.maketrans({c:None for c in ";/"})
    cleaned_block = comb_block.translate(clean).split("}{")
    
    # get longest block
    get_length = lambda x: len(x.split(","))
    return max(map(get_length, cleaned_block))

def check_flatten(x):
    """
    Flatten if it's a list of lists or tuples.
    """
    try:
        x = list(flatten(*x))
    except:
        pass

    return x

def vector2upper_tri_matrix(v):
    # get dimensions of matrix
    n = int((1 + math.sqrt(8 * len(v) + 1)) / 2)
    m = np.zeros((n, n))

    # fill upper triangle
    m[np.triu_indices(n, k=1)] = v

    return m

def vector2list(v):
    """
    Converts the vectorized form of a distance matrix (in the form of 
    a numpy array) to a flattened distance matrix (in the form of a
    list) filled with Decimals. 
    """

    m = vector2upper_tri_matrix(v)

    # fill lower triangle, convert to list, convert elements to Decimal
    return list(map(lambda x: to_decimal(x), (m + m.T).flatten()))

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

def mod_root(root, tree):
    """
    Converts tree into rooted/unrooted tree. Currently being implemented
    in DendroPy. 
    """
    tree.is_rooted = root
    tree.update_bipartitions()
    return tree

def reconstruct_trees(taxa, vectors, true_file):
    """
    Creates trees using DendroPy's Neighbor Joining and UPGMA methods.
    """

    def get_pdm(vector):
        """
        Creates a PhylogeneticDistanceMatrix from the vectorized
        distance matrix.
        """
        csv_string = ",{}\n".format(",".join(taxa.labels()))
        
        m = vector2upper_tri_matrix(vector)

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

    # make distance matrices
    dist_matrices = list(map(get_pdm, vectors))
    nj = lambda x: x.nj_tree()
    upgma = lambda x: x.upgma_tree()

    # make heatmaps
    data_folder, true_name = os.path.split(true_file)
    batch_folder, _ = os.path.split(data_folder)
    dropped_name = true_name[:-9] + ".txt"
    exp_folder = os.path.split(batch_folder)[0]
    heatpath = os.path.join(exp_folder, 'heatmaps')

    if (true_name[1:7] == "200000" and # c ratio is 20
        batch_folder[-1] == "2"):      # number of species is 2
        for i, vect in enumerate(vectors):
            if random.random() < 0.01: # roughly 1/100 become heatmaps
                
                # convert vector to distance matrix
                upper = vector2upper_tri_matrix(vect)
                dist = upper + upper.T

                # write distance matrix to heatmaps folder for later plotting
                plotname = "{}_{}__{}.gz".format(dropped_name[:-ext_len], 
                                              i,
                                              "_".join(taxa.labels()))
                plotpath = os.path.join(heatpath, plotname)
                np.savetxt(plotpath, dist)

    # manually root/unroot
    nj_trees = list(map(partial(mod_root, False), map(nj, dist_matrices)))
    upgma_trees = list(map(partial(mod_root, True), map(upgma, dist_matrices)))

    return nj_trees, upgma_trees

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
        index = to_decimal(index)

        return x[ind] + index%1 * (x[ind+1] - x[ind]) if index % 1 else x[ind]

    return percentile

def param_value(s, tag):
    """
    Finds the value of the parameter specified by 'tag' in 's'.
    """
    params = os.path.basename(s).split("_")
    correct_params = list(filter(lambda tags: tag in tags, params))
    
    assert len(correct_params) == 1
    
    return correct_params[0][len(tag):]

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
        species_depth = int(param_value(sol_file, 'd'))
        theoretical_max = 2 * species_depth 
        nrmse = rmse / theoretical_max

        # calculate percentiles
        p = percentiles_of(imp_err)

        retval = [p(0), p(25), p(50), p(75), p(100), rmse, nrmse]

    elif mode == 'tree': 
        retval = rmse

    return retval

def leaf_stats(solution_file, true_file, batch_folder, dists=[]):
    """
    Finds quartiles, Root Mean Squared Error (RMSE), and Normalized
    Root Mean Squared Error (NRMSE) using files with solution vectors
    and actual distance vectors.
    """
    # read file
    logger = logging.getLogger("leaf_stats")
    try:
        sol_list = np.loadtxt(solution_file)
        true_list = np.loadtxt(true_file)

        assert sol_list.size != 0

    except AssertionError as e:
        logger.warning("Solution file was empty. Please check for " +
                       "imputation errors.")
        raise e

    except FileNotFoundError as e:
        if not sol_list:
            ftype = "Solution File"
            fname = solution_file
        else:
            ftype = "True File"
            fname = true_file
        logger.warning("{} file {} is missing".format(ftype, fname))
        return []

    except OSError as e:
        logger.critical("Could not open either " +
                        "'{}' or '{}'".format(solution_file, true_file) + 
                        ". Please double check the integrity of these files.")
        return []

    # get lists of distances
    imp_dists = list(iterflatten(map(vector2list, sol_list)))
    og_dists = list(iterflatten(map(vector2list, true_list)))

    # calculations
    err = list(map(lambda x: x[0] - x[1], zip(imp_dists, og_dists)))

    return calc(err, sol_file=solution_file), err

def tree_stats(solution_file, true_file, tree_gen, batch_folder, dists):  
    """
    Compiles the BHV and Robinson-Foulds distances, the codimension of
    the BHV distance, and the proportion of well-separated trees.
    """

    # read files
    logger = logging.getLogger("tree_stats")
    try:
        sol_list = np.loadtxt(solution_file)
        true_list = np.loadtxt(true_file)
        sep_file = os.path.join(batch_folder, "nexus/separated.txt")
        og_sep = [float(np.loadtxt(sep_file))]

        assert sol_list.size != 0

    except AssertionError as e:
        logger.warning("Solution file was empty. Please check for " +
                       "imputation errors.")
        raise e

    except FileNotFoundError as e:
        if not sol_list:
            ftype = "Solution File"
            fname = solution_file
        elif not true_list:
            ftype = "True File"
            fname = true_file
        else:
            ftype = "Well-Separation"
            fname = sep_file

        logger.warning("{} file {} is missing".format(ftype, fname))

        return []

    # get trees
    imp_nj, imp_upgma = tree_gen(sol_list, true_file)
    og_nj, og_upgma = tree_gen(true_list, true_file)

    nj_trees = list(zip(imp_nj, og_nj))
    upgma_trees = list(zip(imp_upgma, og_upgma))
    sep_trees = [imp_nj, imp_upgma, og_nj, og_upgma]

    # let bhv_r and bhv_u be mapped
    bhv_r = partial(bhv_dist, batch_folder=batch_folder, rooted=True)
    bhv_u = partial(bhv_dist, batch_folder=batch_folder, rooted=False)

    # distance calculations
    tree_dists = []
    codim = []
    if 'bhv' in dists:
        bhv_nj, codim_nj = list(zip(*list(map(bhv_u, nj_trees))))
        bhv_upgma, codim_upgma = list(zip(*list(map(bhv_r, upgma_trees))))
        tree_dists += [bhv_nj, bhv_upgma]
        codim = [codim_nj, codim_upgma]
    if 'rf' in dists:
        rf_nj = list(map(rf_dist, nj_trees))
        rf_upgma = list(map(rf_dist, upgma_trees))
        tree_dists += [rf_nj, rf_upgma]
    if 'norm' in dists:
        l2_nj = list(map(l2_norm, nj_trees))
        l2_upgma = list(map(l2_norm, upgma_trees))
        tree_dists += [l2_nj, l2_upgma]

    # proportion of good clades
    well_sep = list(map(prop_separated_trees, sep_trees)) + og_sep
    
    # rearrange 
    tree_dists_plus = tree_dists + codim

    # make sure all fields are the same length
    assert all(map(lambda x: len(x) == len(nj_trees), tree_dists_plus))

    # make summaries and return 
    tree_calc = partial(calc, mode="tree")
    return list(map(tree_calc, tree_dists)) + well_sep, tree_dists_plus

def write_stats(desc, get_stats, tup, dists=[]):
    """
    Writes data to the relevant file in the stats subdirectory.

    'get_stats' is a function which takes as input an exploded 'tup' and
    returns the stats to write to file. 
    'desc' is a description tag for the filenames.
    """

    solution_file = tup[0]
    batch_folder = tup[-1]
    data = get_stats(*tup, dists)

    assert len(data) == 2

    stats = os.path.join(batch_folder, "stats")
    if not os.path.isdir(stats):
        os.makedirs(stats)

    basename = os.path.basename(solution_file)[:-ext_len] + desc

    # separate and flatten data
    summary, all_data = list(map(check_flatten, data))

    # write data
    with gzip.open(os.path.join(stats, basename + '.txt.gz'), 'wt') as f:
        f.write(" ".join(map(str, summary)))
    with gzip.open(os.path.join(stats, basename + "_all.txt.gz"), 'wt') as f:
        f.write(" ".join(map(str, all_data)))

def match(infile, file_list, tags):
    """
    Finds files in file_list containing the same tag value(s) as infile.
    """
    if type(tags) is list:
        tag_values = list(filter(lambda x: any(tag in x for tag in tags), 
                                 infile[:-ext_len].split("_")))
    else:
        tag_values = list(filter(lambda x: tags in x, infile[:-ext_len].split("_")))
    
    match_files = list(filter(lambda f: all(tag in f for tag in tag_values), 
                              file_list))

    try:
        assert len(match_files) == 1
    except AssertionError as e:
        logger = logging.getLogger("match")
        logger.debug(infile)
        logger.debug(file_list)
        logger.debug(tags)
        logger.debugs(match_files)
        logger.exception('')
        raise e

    return match_files[0]

def analyze(batch_folder, dists):
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
    true_files = tagged_files("data", "_true.txt.gz")
    sol_files = tagged_files("solutions", ".sol.gz")
    nexus_files = tagged_files("nexus", "_e")

    # match solution files to files from which they were created
    match_sol_true = lambda sol: match(sol, true_files, ['e', 'p'])

    # match solution files with a tree generation function
    sol2nexus = {sol: match(sol, nexus_files, 'e') for sol in sol_files}
    nex2trees = {n: partial(reconstruct_trees, make_taxa(n)) for n in nexus_files}
    match_sol_treegen = lambda sol: nex2trees[sol2nexus[sol]]

    # make functions to create and write stats
    write_leaves = partial(write_stats, '', leaf_stats)
    write_trees = partial(write_stats, '_tree', tree_stats, dists=dists)

    # make tuples to create and write stats 
    stats_files = zip(sol_files, 
                      map(match_sol_true, sol_files),
                      [batch_folder]*len(sol_files))
    tree_files = zip(sol_files,
                     map(match_sol_true, sol_files), 
                     map(match_sol_treegen, sol_files),
                     [batch_folder]*len(sol_files))

    # create and write stats
    list(map(write_leaves, stats_files))
    list(map(write_trees, tree_files))

def values_from_tags(tags):
    """
    Returns a list of lists of the values of each parameter in tags.
    """

    # check if string is a valid parameter
    valid_param = lambda s: s[0].isalpha() and s[1].isdigit()
    # split by underscore
    def split(s):
        if any(map(lambda x: x.isdigit(), s[-ext_len:])):
            return list(filter(lambda x: x, s.split("_")))
        else:
            return list(filter(lambda x: x, s[:-ext_len].split("_")))

    param_list = set(filter(valid_param, iterflatten(map(split, tags))))  

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

def get_row(vect, param_values, exp_folder, rowsize, dists, modifier=""):
    """
    Get the row of values to store in the CSV file according to
    vect, a vector describing the values of parameters in the
    experiment.
    """

    # get values of current parameters
    get_vals = lambda x: param_values[x[0]][x[1]]
    params = list(map(get_vals, enumerate(vect)))

    # get paths
    batch_folder = batch_analyze.format(*params)
    nameroot = fileroot.format(*params)
    filename = nameroot + modifier + ext
    stats_path = os.path.join(exp_folder,
                              batch_folder,
                              "stats",
                              filename)

    # get data
    logger = logging.getLogger("get_row")
    try:
        # stats_path points to non-tree_all stats file
        data = list(np.loadtxt(stats_path))

        row = params + data 

        if not data:
            row = [np.nan] * rowsize
        # data is for tree_all and must be padded depending on the
        # number of trees 
        elif rowsize != len(row):
            # split data into num_pieces=len(ex_tree_types) pieces
            num = {'bhv': 4, 'rf': 2, 'norm': 2}
            nzr = {'bhv': 2, 'rf': 2, 'norm': 0}
            num_pieces = sum(num[dist] for dist in dists)
            k = int(len(data)/num_pieces)
            nzr = sum(nzr[dist] for dist in dists)

            split_k = [data[i:i+k] for i in range(0, len(data), k)]

            # pad each piece
            pad_len = int((rowsize - len(params) - len(data) - nzr)/num_pieces)
            padded = [q + [float('nan')] * pad_len for q in split_k]

            # find zero ratio
            zero = [len([x for x in q if x == 0]) for q in split_k[:nzr]]
            zr = [z/len(total) for z, total in zip(zero, split_k[:nzr])]
            
            # create row
            row = params + list(flatten(*padded)) + zr

            assert len(row) == rowsize

    except ValueError as e:
        logger.debug(stats_path)
        logger.debug("Filename: {}".format(filename))
        logger.debug("k: {}".format(k))
        raise e

    except AssertionError as e:
        logger.debug("Filename: {}".format(filename))
        logger.debug("rowsize: {}, len(row): {}".format(rowsize, len(row)))
        logger.debug("pad_len: {}, k: {}".format(pad_len, k))
        logger.debug("data: {}".format(list(map(len, split_k))))
        raise e

    except FileNotFoundError as e:
        logger.warning("Stats file {} not found. ".format(stats_path) +
                       "If the corresponding solution file exists, " +
                       "please run the 'analyze' step again. " +
                       "Otherwise, double check the imputation step. " +
                       "This warning indicates possibly biased " +
                       "missingness in your data.")
        row = [np.nan] * rowsize

    return row

def compile_stats(exp_folder, dists):
    """
    Creates a csv file with summary statistics from each batch_folder in 
    exp_folder.
    """

    # make sure exp_folder is absolute path
    exp_folder = os.path.abspath(exp_folder)

    # figure out which experiments were run
    valid = lambda name: name[0].isalpha() and name[1].isdigit()
    supertags = list(filter(valid, os.listdir(exp_folder)))
    subtags = list(filter(valid, os.listdir(os.path.join(exp_folder, 
                                                         supertags[0],
                                                         'stats'))))

    # get values of parameters used in this experiment
    param_values = values_from_tags(supertags + subtags)

    # get number of gene trees
    n_gene_trees = max(map(lambda x: int(x[1:]),
                                filter(lambda x: 'g' in x, 
                                       iterflatten(map(lambda x: x.split("_"),
                                                       supertags)))))

    # helper functions
    expand_cols = lambda x: [x]*n_gene_trees

    # specify column names
    base_cols = ["Species Depth", "Species Tree Index", "Number of Gene Trees", 
                 "Number of Individuals per Species", "Method",
                 "Effective Population Size", "Leaf Dropping Probability",
                 "Number of Species"]
    stats_types = ["Abs Imputation Error (AIE) min", 
                   "AIE lower quartile", "AIE median", "AIE upper quartile",
                   "AIE max", "Imputation RMSE", "Imputation NRMSE"]

    base_tree_types = []
    codims = []
    zr = []
    if 'bhv' in dists:
        base_tree_types += ['bhv_nj', 'bhv_upgma']
        zr += ["zr_bhv_nj", "zr_bhv_upgma"]
        codims = ["nj_codim", "upgma_codim"]
    if 'rf' in dists:
        base_tree_types += ['rf_nj', 'rf_upgma']
        zr += ["zr_rf_nj", "zr_rf_upgma"]
    if 'norm' in dists:
        base_tree_types += ['l2_nj', 'l2_upgma']

    tree_types = base_tree_types + ["Proportion Separated Imputed NJ Trees",
                                    "Proportion Separated Imputed UPGMA Trees",
                                    "Proportion Separated Reconstructed NJ Trees",
                                    "Proportion Separated Reconstructed UPGMA Trees",
                                    "Proportion Separated Original Trees"]
    ex_tree_types = base_tree_types + codims

    expanded_tree_types = list(iterflatten(map(expand_cols, ex_tree_types)))
    tree_all_types = expanded_tree_types + zr

    # define columns
    extras = [stats_types, tree_types, tree_all_types]
    cols = [base_cols + extra for extra in extras]
    stats_cols, tree_cols, tree_all_cols = cols 

    # get sizes to create the numpy array
    dimensions = list(map(len, param_values))
    size = reduce(operator.mul, dimensions, 1)

    # build data arrays
    build_array = lambda cols: np.empty((size, len(cols)), dtype=float)
    data, tree, tree_all = list(map(build_array, cols))

    # fill data arrays
    coordinates = itertools.product(*list(map(range, dimensions)))
    for ind, coordinate in enumerate(coordinates):
        data[ind] = get_row(coordinate,
                            param_values,
                            exp_folder,
                            len(cols[0]),
                            dists, 
                            "")
        tree[ind] = get_row(coordinate, 
                            param_values, 
                            exp_folder,
                            len(cols[1]),
                            dists,
                            "_tree")
        tree_all[ind] = get_row(coordinate, 
                                param_values,
                                exp_folder,
                                len(cols[2]),
                                dists,
                                "_tree_all")

    # make dataframes
    stats = pd.DataFrame(data, columns=stats_cols)
    tree = pd.DataFrame(tree, columns=tree_cols)
    tree_all = pd.DataFrame(tree_all, columns=tree_all_cols)
    csvs = [stats, tree, tree_all]

    # merge old csv files
    cnames = ["interleaf_error.csv", 
              "intertree_error.csv",
              "intertree_all.csv"]
    fullnames = list(map(lambda x: os.path.join(exp_folder, x),
                         cnames))
    existing = [csv in os.listdir(exp_folder) for csv in cnames]
    for i, csv in enumerate(fullnames):
        if not existing[i]: continue
        old = pd.read_csv(csv)
        old_cols = list(old.columns.values)
        new_cols = list(csvs[i].columns.values)
        dif_cols = ['' if c in new_cols else c for c in old_cols]
        dif_cols = list(filter(lambda x: bool(x), dif_cols))
        csvs[i] = pd.concat([csvs[i], old[dif_cols]], axis=1)

    # write csvs
    write_csv = lambda tup: tup[0].to_csv(tup[1])
    list(map(write_csv, zip(csvs, fullnames)))

def make_plots(exp_folder):
    # write heatmaps
    sns.set_palette('colorblind')
    exp_folder = os.path.realpath(exp_folder)
    heatpath = os.path.join(exp_folder, 'heatmaps')
    fig = plt.figure()

    for plotname in os.listdir(heatpath):
        if plotname[-4:] == '.png' or plotname[0] == '.': continue
        plotpath = os.path.join(heatpath, plotname)
        dist = np.loadtxt(plotpath)

        name, labelstring = plotpath.split('__')
        labels = labelstring.split("_")

        # file gymnastics
        outname = os.path.join(heatpath, name)
        ind = name.split("_")[-1]
        _, plot_name = os.path.split(name[:-(len(ind)+1)])
        valid = lambda name: name[0].isalpha() and name[1].isdigit()
        supertags = list(filter(valid, os.listdir(exp_folder)))
        b_len = len(supertags[0].split('_'))
        batch_base = "_".join(os.path.split(plot_name)[1].split('_')[:b_len])
        data_name = os.path.join(exp_folder,
                                 batch_base, 
                                 "data",
                                 plot_name + ext)

        # extract distance matrix for tree with dropped leaves
        with gzip.open(data_name, 'rt') as f:
            lines = (line for line in f if len(line.split()) > 2)
            all_dropped = np.loadtxt(lines)
            dropped = vector2upper_tri_matrix(all_dropped[int(ind)])
            dropped += dropped.T

        # convert -1 to -max in the dropped matrix for visualization
        dmax = max(np.amax(dropped), np.amax(dist))
        dropped = np.vectorize(lambda x: -dmax if x == -1 else x)(dropped)

        # matrix of tree with dropped leaves
        dhm = pd.DataFrame(dropped, columns=labels, index=labels)
        dhm = dhm.reindex_axis(sorted(dhm.columns), axis='columns')
        dhm = dhm.reindex_axis(sorted(dhm.index), axis='index')

        # matrix of tree with imputed leaves
        ihm = pd.DataFrame(dist, columns=labels, index=labels)
        ihm = ihm.reindex_axis(sorted(ihm.columns), axis='columns')
        ihm = ihm.reindex_axis(sorted(ihm.index), axis='index')

        # dropped heatmap
        sns.heatmap(dhm)
        fig.savefig(outname + "_dropped")
        fig.clear()

        # imputed heatmap
        sns.heatmap(ihm, vmin=-dmax)
        fig.savefig(outname + "_imputed")
        fig.clear()

        # remove heatmap data file
        os.remove(plotpath)
        
    plt.close()

    # run R plots for descriptive statistics
    summary_graphs = "Rscript_$_summary.R_$_{}".format(exp_folder)
    summary = get_output(summary_graphs.split("_$_"), ignore_error=True)
    outlier_counts = "Rscript_$_outliers.R_$_{}".format(exp_folder)
    outlier = get_output(outlier_counts.split("_$_"), ignore_error=True)
    if summary == "error":
        logging.warning("\nError while making summary graphs. This " +
                        "error does not affect the output csvs.")
    if outlier == "error":
        logging.warning("\nError while making outlier graphs. This " +
                        "error does not affect the output csvs.")
        try:
            os.remove("Rplots.pdf")
        except FileNotFoundError as e:
            pass

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    exp_folder = args['<exp_folder>']
    dists = args['<dists>']

    compile_stats(exp_folder, dists)