"""
WARNING: Do not edit any files that this program relies on while it
is running. You will end up changing the experimental environment!

To safely run, copy this file and other files into another folder and 
run from that other location. Then you may edit the original files.
"""

from subprocess import check_call as cc
import itertools

def setup(name_formula, names):

    def one_run(tup):
        cc(name_formula.format(**dict(zip(names, tup))), shell=True)

    return one_run

if __name__ == "__main__":

    # edit these lists to your heart's desire
    c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20].reverse()
    genes = [10, 20, 30, 40, 50, 60]
    inds = [4, 5, 6, 10]
    methods = [1, 2]
    probs = [0.1, 0.2, 0.5]
    species = [2, 3, 4, 5]

    # put together the python call
    root = "python3 pipeline.py generate {}"
    opts = " {c} -g{genes} -i{inds} -m{method} -p{prob} -s{sp}"
    
    # names of the fields used in the format string, in the same order
    # as they appear in the itertools.product tuple.
    names = ["c", "genes", "inds", "method", "prob", "sp"]
    
    # feel free also to change the experiment name path
    name = "../experiments/c{c}_g{genes}_i{inds}_m{method}_p{prob}_s{sp}"

    list(map(setup(root.format(name + opts), names),
             itertools.product(c, genes, inds, methods, probs, species)))