"""
Usage: analyze.py <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders.
"""

import decimal
import docopt
from itertools import chain
import numpy as np
import math
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
    imputed = list(chain.from_iterable(map(vector2matrix, vec_list)))
    original = list(chain.from_iterable(map(vector2matrix, true_list)))

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

    def write_stats(x):
        with open(os.path.join(stats, basename+'_stats.txt'), 'w') as f:
            list(map(lambda i: f.write(str(i)+" "), x)) 
    return write_stats

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

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    infolder = args['<infolder>']

    anlyze(infolder)