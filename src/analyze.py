"""
Usage: analyze.py <infolder>

Options:
  -h, --help            Show this help message.
  <infolder>            Folder containing data, nexus, and solutions
                        subfolders.
"""

import docopt
import numpy as np
import math
import os

def vector2matrix(v):
    n = int((1 + math.sqrt(8 * len(v) + 1)) / 2)
    m = np.zeros((n, n))
    m[np.triu_indices(n, k=1)] = v
    return m + m.T

def get_stats(x):
    solutions, true = x
    # read file
    vec_list = np.loadtxt(solutions)
    true_list = np.loadtxt(true)

    if vec_list.size == 0:
        return ["Infile was empty. Please check for imputation errors."]

    # make distace matrix array
    imputed = np.array(list(map(vector2matrix, vec_list)))
    original = np.array(list(map(vector2matrix, true_list)))

    sq_err = (imputed - original)**2
    sse = np.sum(sq_err)/np.count_nonzero(sq_err)
    mse = sse/imputed.size
    rmse = math.sqrt(mse)

    return sse, mse, rmse

def write_stats(infolder):
    def write(x):
        with open(os.path.join(infolder, "stats.txt"), 'w') as f:
            list(map(lambda i: f.write(str(i)+"\t"), x)) 
    return write

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    infolder = args['<infolder>']

    inf_path = lambda x: os.path.join(infolder, x)
    files = lambda x: map(lambda y: os.path.join(infolder, x, y),
                          os.listdir(inf_path(x)))

    true_files = filter(lambda x: "_true.txt" in os.path.basename(x), 
                        files("data"))
    solution_files = filter(lambda x: ".sol" in os.path.basename(x),
                            files("solutions"))

    analyze = lambda x: write_stats(infolder)(get_stats(x))

    list(map(analyze, zip(solution_files, true_files)))