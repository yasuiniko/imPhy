"""
Usage: analyze.py <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders.
"""

import decimal
import docopt
import itertools
import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import os

def vector2matrix(v):
    n = int((1 + math.sqrt(8 * len(v) + 1)) / 2)
    m = np.zeros((n, n))
    m[np.triu_indices(n, k=1)] = v
    return list(map(lambda x: decimal.Decimal(str(x)), (m + m.T).flatten()))

def percentiles_of(x):
    x = sorted(x)
    n = len(x)
    
    def percentile(p):
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

def get_stats(x):
    solutions, true = x
    # read file
    vec_list = np.loadtxt(solutions)
    true_list = np.loadtxt(true)

    if vec_list.size == 0:
        return ["Infile was empty. Please check for imputation errors."]

    to_decimal = lambda x: decimal.Decimal(str(x))

    # make distace matrix array
    imputed = list(itertools.chain.from_iterable(map(vector2matrix, vec_list)))
    original = list(itertools.chain.from_iterable(map(vector2matrix, true_list)))

    sq_err = list(map(lambda x: (x[0] - x[1])**2, zip(imputed, original)))
    imp_sq_err = list(filter(lambda x: x != 0, sq_err))
    imp_err = list(map(lambda x: x.sqrt(), imp_sq_err))
    sse = to_decimal(sum(imp_sq_err))
    mse = sse/len(imp_sq_err)
    rmse = math.sqrt(mse)

    p = percentiles_of(imp_err) 

    return p(0), p(25), p(50), p(75), p(100), rmse

def write_to(infolder, sol):

    stats = os.path.join(infolder, "stats")
    if not os.path.isdir(stats):
        os.makedirs(stats)

    basename = os.path.basename(sol)[:-4]

    def write(x):
        with open(os.path.join(stats, basename+'_stats.txt'), 'w') as f:
            list(map(lambda i: f.write(str(i)+" "), x)) 
    return write

def analyze(infolder):
    inf_path = lambda x: os.path.join(infolder, x)
    files = lambda x: map(lambda y: os.path.join(infolder, x, y),
                          os.listdir(inf_path(x)))

    true_files = filter(lambda x: "_true.txt" in os.path.basename(x), 
                        files("data"))
    solution_files = filter(lambda x: ".sol" in os.path.basename(x),
                            files("solutions"))

    write_stats = lambda x: write_to(infolder, x[0])(get_stats(x))

    list(map(write_stats, zip(solution_files, true_files)))

def summary(expfolder):
    from_iterable = itertools.chain.from_iterable
    expfolder = os.path.abspath(expfolder)

    # figure out which experiments were run
    letter_number = lambda x: x[0].isalpha() and x[1].isdigit()
    split = lambda x: x.split("_")
    supertrials = list(filter(letter_number, os.listdir(expfolder)))
    subtrials = filter(letter_number, os.listdir(os.path.join(expfolder, 
                                                 supertrials[0],
                                                 'stats')))
    superlabels = set(filter(letter_number,
                             from_iterable(map(split, supertrials))))
    sublabels = set(filter(letter_number, 
                           from_iterable(map(split, subtrials))))
    labels = superlabels | sublabels
    vals = lambda x: list(map(lambda label: label.split(x)[-1],
                              filter(lambda label: x in label, labels)))
    c = vals("c")
    g = vals("g")
    i = vals("i")
    m = vals("m")
    p = vals("p")
    s = vals("s")
    dimensions = list(map(len, (c, g, i, m, p, s)))
    size = len(c) * len(g) * len(i) * len(m) * len(p) * len(s)

    def get_row(args):
        infolder, filename, c = args
        nums = filter(lambda x: x not in "gimps", os.path.basename(infolder))
        info = list(map(float, ''.join(nums).split("_")))
        filepath = os.path.join(infolder, "stats", filename)

        try: 
            row = [c] + info + list(map(float, np.loadtxt(filepath)))
        except ValueError:
            row = [float('nan')]*12

        return row

    def get_folder(tup):
        folder_name = "g{}_i{}_m{}_p{}_s{}".format(*tup[1:])
        filename = "g{1}_i{2}_m{3}_p{4}_s{5}_c{0}_genes_1_stats.txt".format(*tup)
        return os.path.join(expfolder, folder_name), filename, tup[0]

    def get_labels(tup):
        return c[tup[0]], g[tup[1]], i[tup[2]], m[tup[3]], p[tup[4]], s[tup[5]]

    # build data array
    data = np.empty((size, 12), dtype=float)
    for ind, index in enumerate(itertools.product(*list(map(range, dimensions)))):
        data[ind] = get_row(get_folder(get_labels(index)))
    
    # put data in a pandas DataFrame
    cols = ["SD/Ne", "Number of Gene Trees", 
            "Number of Individuals per Species", "Method", 
            "Leaf Dropping Probability", "Number of Species", 
            "IE min", "IE lower quartile", "IE median", 
            "IE upper quartile", "IE max", "Imputation RMSE"]
    data = pd.DataFrame(data, columns=cols)
    data.to_csv(os.path.join(expfolder, "summary.csv"))

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    infolder = args['<infolder>']

    summary(infolder)