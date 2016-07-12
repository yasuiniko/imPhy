from subprocess import check_call as cc
import itertools

def one_run(tup, name_formula):
    names = ["c", "genes", "inds", "method", "prob", "sp"]
    cc(name_formula.format(**dict(zip(names, tup))), shell=True)

if __name__ == "__main__":
    c = [0.6, 0.7, 0.8, 0.9, 1, 2, 4, 6, 8, 10, 20]
    genes = [10, 20, 30, 40, 50, 60]
    inds = [4, 5, 6, 10]
    methods = range(1, 2)
    probs = [0.1, 0.2, 0.5]
    species = [2, 3, 4, 5]
    root = "python3 pipeline.py generate {}"
    opts = " {c} -g{genes} -i{inds} -m{method} -p{prob} -s{sp}"
    name = "../experiments/c{c}_g{genes}_i{inds}_m{method}_p{prob}_s{sp}"

    verbose_run = lambda tup: one_run(tup, root.format(name + opts))

    list(map(verbose_run, 
             itertools.product(c, genes, inds, methods, probs, species)))