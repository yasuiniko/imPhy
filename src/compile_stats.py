"""
Usage: compile_stats.py <exp_folder>

Options:
  -h, --help            Show this help message.
  <exp_folder>          Folder containing data, nexus, and solutions
                        subfolders.
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

def rf_dist(trees):
    """
    Calculates the Robinson-Foulds distance between two trees.
    """
    assert len(trees) == 2
    dist = dendropy.calculate.treecompare.symmetric_difference(*trees)
    n = len(trees[0].leaf_nodes())

    return dist/(2*(n-3))

def bhv_dist(trees, batch_folder, rooted):
    """
    Wrapper for Owen and Provan's GTP code to calculate the BHV 
    distance between two trees.
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
        try:
            output = str(subprocess.check_output(gtp.split("^^")), "utf-8")
        except subprocess.CalledProcessError as e:
            print(output)
            raise e

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
    find longest subblock in 'combinatorial types' block
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
    if x and type(x[0]) == list or type(x[0]) == tuple:
        x = list(flatten(*x))
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
    Converts tree into rooted/unrooted tree.
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

    # make trees
    dist_matrices = list(map(get_pdm, vectors))
    nj = lambda x: x.nj_tree()
    upgma = lambda x: x.upgma_tree()

    # make heatmaps
    data_folder, true_name = os.path.split(true_file)
    batch_folder, _ = os.path.split(data_folder)
    dropped_name = true_name[:-9] + ".txt"
    exp_folder = os.path.split(batch_folder)[0]
    heatpath = os.path.join(exp_folder, 'heatmaps')

    if (true_name[1:7] == "200000" and 
        batch_folder[-1] == "2"):
        for i, vect in enumerate(vectors):
            if random.random() < 0.01:
                upper = vector2upper_tri_matrix(vect)
                dist = upper + upper.T
                plotname = "{}_{}__{}".format(dropped_name[:-4], 
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

def leaf_stats(solution_file, true_file, batch_folder):
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
    imp_dists = list(iterflatten(map(vector2list, sol_list)))
    og_dists = list(iterflatten(map(vector2list, true_list)))

    # calculations
    err = list(map(lambda x: x[0] - x[1], zip(imp_dists, og_dists)))

    return calc(err, sol_file=solution_file), err

def tree_stats(solution_file, true_file, tree_gen, batch_folder):  
    """
    Compiles the BHV and Robinson-Foulds distances, the codimension of
    the BHV distance, and the proportion of well-separated trees.
    """

    # read files
    sol_list = np.loadtxt(solution_file)
    true_list = np.loadtxt(true_file)
    sep_file = os.path.join(batch_folder, "nexus/separated.txt")
    og_sep = [float(np.loadtxt(sep_file))]

    # get trees
    imp_nj, imp_upgma = tree_gen(sol_list, true_file)
    og_nj, og_upgma = tree_gen(true_list, true_file)

    nj_trees = list(zip(imp_nj, og_nj))
    upgma_trees = list(zip(imp_upgma, og_upgma))
    sep_trees = [imp_nj, imp_upgma, og_nj, og_upgma]

    # turn bhv_dist into a monad
    bhv_r = partial(bhv_dist, batch_folder=batch_folder, rooted=True)
    bhv_u = partial(bhv_dist, batch_folder=batch_folder, rooted=False)

    # distance calculations
    bhv_nj, codim_nj = list(zip(*list(map(bhv_u, nj_trees))))
    bhv_upgma, codim_upgma = list(zip(*list(map(bhv_r, upgma_trees))))
    rf_nj = list(map(rf_dist, nj_trees))
    rf_upgma = list(map(rf_dist, upgma_trees))

    # proportion of good clades
    well_sep = list(map(prop_separated_trees, sep_trees)) + og_sep
    
    # rearrange 
    tree_dists = [bhv_nj, bhv_upgma, rf_nj, rf_upgma]
    tree_dists_plus = [*tree_dists, codim_nj, codim_upgma]

    # make sure all fields are the same length
    assert all(map(lambda x: len(x) == len(nj_trees), tree_dists_plus))

    # make summaries and return 
    tree_calc = partial(calc, mode="tree")
    return list(map(tree_calc, tree_dists)) + well_sep, tree_dists_plus

def write_stats(desc, fun, tup):
    """
    Writes data to the relevant file in the stats subdirectory.

    'fun' is a data generating function which takes as input an 
    exploded 'tup'. desc is a description tag for the filenames.
    """

    solution_file = tup[0]
    batch_folder = tup[-1]
    data = fun(*tup)

    stats = os.path.join(batch_folder, "stats")
    if not os.path.isdir(stats):
        os.makedirs(stats)

    basename = os.path.basename(solution_file)[:-4] + desc

    assert len(data) == 2

    # separate and flatten data
    summary, all_data = list(map(check_flatten, data))
    
    # write data
    with open(os.path.join(stats, basename + '.txt'), 'w') as f:
        f.write(" ".join(map(str, summary)))
    with open(os.path.join(stats, basename + "_all.txt"), 'w')as f:
        f.write(" ".join(map(str, all_data)))

def match(infile, file_list, tags):
    """
    Finds files in file_list containing the same tag value(s) as infile.
    """
    if type(tags) is list:
        tag_values = list(filter(lambda x: any(tag in x for tag in tags), 
                                 infile[:-4].split("_")))
    else:
        tag_values = list(filter(lambda x: tags in x, infile[:-4].split("_")))
    
    match_files = list(filter(lambda f: all(tag in f for tag in tag_values), 
                              file_list))

    try:
        assert len(match_files) == 1
    except AssertionError as e:
        print(infile)
        print(file_list)
        print(tags)
        print(match_files)
        raise e

    return match_files[0]

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
    match_sol_true = lambda sol: match(sol, true_files, ['e', 'p'])

    # match solution files with a tree generation function
    sol2nexus = {sol: match(sol, nexus_files, 'e') for sol in sol_files}
    nex2trees = {n: partial(reconstruct_trees, make_taxa(n)) for n in nexus_files}
    match_sol_treegen = lambda sol: nex2trees[sol2nexus[sol]]

    # make functions to create and write stats
    write_leaves = partial(write_stats, '', leaf_stats)
    write_trees = partial(write_stats, '_tree', tree_stats)

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
        if list(filter(lambda x: x.isdigit(), s[-4:])):
            return list(filter(lambda x: x, s.split("_")))
        else:
            return list(filter(lambda x: x, s[:-4].split("_")))

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

def get_row(vect, param_values, exp_folder, rowsize, modifier=""):
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
    filename = nameroot + modifier + ".txt"
    stats_path = os.path.join(exp_folder,
                              batch_folder,
                              "stats",
                              filename)

    # get data
    try:
        # stats_path points to non-tree stats file
        data = list(np.loadtxt(stats_path))
        row = params + data 

        # data is for tree_all and must be padded depending on the
        # number of trees 
        if rowsize != len(row):
            # split data into num_pieces pieces 
            num_pieces = 6
            k = int(len(data)/num_pieces)

            split_k = [data[i:i+k] for i in range(0, len(data), k)]

            # pad each piece
            pad_len = int((rowsize - len(params) - len(data) - 4)/num_pieces)
            padded = [q + [float('nan')] * pad_len for q in split_k]

            # find zero ratio
            zero = [len([x for x in q if x == 0]) for q in split_k[:4]]
            zr = [z/len(total) for z, total in zip(zero, split_k[:4])]
            
            # create row
            row = params + list(flatten(*padded)) + zr

            assert len(row) == rowsize

    except ValueError as e:
        print("Filename: {}".format(filename))
        print("k: {}".format(k))
        raise e
        
        row = [float('nan')] * rowsize

        # # code for avoiding lists, but should avoid lists when writing
        # with open(stats_path, 'r') as f:
        #     cln = lambda x: x.translate({ord(c): None for c in '[],'})
        #     data = list(map(cln, f.readline().split()))
        # if list(filter(lambda x: x.isdigit(), data)):
        #     row = params + data
        # else:
        #     row = [float('nan')] * rowsize

    except AssertionError as e:
        print("Filename: {}".format(filename))
        print("rowsize: {}, len(row): {}".format(rowsize, len(row)))
        print("pad_len: {}, k: {}".format(pad_len, k))
        print("data: {}".format(list(map(len, split_k))))
        raise e

    except FileNotFoundError as e:
        print("Please re-run the 'analyze' step in experiment.py.")
        raise e

    return row

def compile_stats(exp_folder):
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
    base_tree_types = ["bhv_nj", "rf_nj", "bhv_upgma", "rf_upgma"]
    tree_types = base_tree_types + ["Proportion Separated Imputed NJ Trees",
                                    "Proportion Separated Imputed UPGMA Trees",
                                    "Proportion Separated Reconstructed NJ Trees",
                                    "Proportion Separated Reconstructed UPGMA Trees",
                                    "Proportion Separated Original Trees"]
    ex_tree_types = base_tree_types + ["nj_codim", "upgma_codim"]

    expanded_tree_types = list(iterflatten(map(expand_cols, ex_tree_types)))
    zr = ["zr_bhv_nj","zr_rf_nj",
           "zr_bhv_upgma","zr_rf_upgma",]
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
                            "")
        tree[ind] = get_row(coordinate, 
                            param_values, 
                            exp_folder,
                            len(cols[1]),
                            "_tree")
        tree_all[ind] = get_row(coordinate, 
                                param_values,
                                exp_folder,
                                len(cols[2]),
                                "_tree_all")

    # write stats summary
    outpath = os.path.join(exp_folder, "summary.csv")
    pd.DataFrame(data, columns=stats_cols).to_csv(outpath)
    # write tree summary
    outpath = os.path.join(exp_folder, "tree_summary.csv")
    pd.DataFrame(tree, columns=tree_cols).to_csv(outpath)
    # write tree all
    outpath = os.path.join(exp_folder, "tree_all.csv")
    pd.DataFrame(tree_all, columns=tree_all_cols).to_csv(outpath)

    # write heatmaps
    sns.set_palette('colorblind')
    heatpath = os.path.join(exp_folder, 'heatmaps')
    fig = plt.figure()
    for plotname in os.listdir(heatpath):
        if plotname[-4:] == '.png' or plotname[0] == '.': continue
        plotpath = os.path.join(heatpath, plotname)
        dist = np.loadtxt(plotpath)

        name, labelstring = plotpath.split('__')
        labels = labelstring.split("_")

        outname = os.path.join(heatpath, name)
        ind = name.split("_")[-1]
        _, plot_name = os.path.split(name[:-(len(ind)+1)])
        b_len = len(supertags[0].split('_'))
        batch_base = "_".join(os.path.split(plot_name)[1].split('_')[:b_len])
        data_name = os.path.join(exp_folder,
                                 batch_base, 
                                 "data",
                                 plot_name + ".txt")
        with open(data_name, 'r') as f:
            lines = (line for line in f if len(line.split()) > 2)
            all_dropped = np.loadtxt(lines)
            dropped = vector2upper_tri_matrix(all_dropped[int(ind)])
            dropped += dropped.T

        dmax = max(np.amax(dropped), np.amax(dist))
        set_min = np.vectorize(lambda x: -dmax if x == -1 else x)
        dropped = set_min(dropped)

        dhm = pd.DataFrame(dropped, columns=labels, index=labels)
        dhm = dhm.reindex_axis(sorted(dhm.columns), axis='columns')
        dhm = dhm.reindex_axis(sorted(dhm.index), axis='index')

        ihm = pd.DataFrame(dist, columns=labels, index=labels)
        ihm = ihm.reindex_axis(sorted(ihm.columns), axis='columns')
        ihm = ihm.reindex_axis(sorted(ihm.index), axis='index')

        sns.heatmap(dhm, yticklabels=labels)
        fig.savefig(outname + "_dropped")
        fig.clear()

        sns.heatmap(ihm, yticklabels=labels, vmin=-dmax)
        fig.savefig(outname + "_imputed")
        fig.clear()

        os.remove(plotpath)
        
    plt.close()

    make_graphs = "Rscript_$_summary.R_$_{}".format(exp_folder)
    subprocess.check_call(make_graphs.split("_$_"))

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    exp_folder = args['<exp_folder>']

    compile_stats(exp_folder)