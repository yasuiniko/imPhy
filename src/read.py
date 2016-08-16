import numpy as np
import scipy
import scipy.spatial
import math
import itertools
import functools

def vector2matrix(v):
    n = int((1 + math.sqrt(8 * len(v) + 1)) / 2)
    
    m = np.zeros((n, n))

    m[np.triu_indices(n, k=1)] = v

    return m + m.T

def compose(*functions):
    return functools.reduce(lambda f, g: lambda x: f(g(x)), 
                            functions,
                            lambda x: x)

if __name__ == '__main__':

    # quality of life
    np.set_printoptions(precision=2, suppress=True, linewidth=110)

    # setup
    filename = '/Users/nikoyasui/Desktop/Japan/trees/phyclust/src/test.txt'

    # read file
    with open(filename, 'r') as f:
        dlist = []
        n_sp, n_leaves = list(map(int, f.readline().strip().split(' ')))
        n_inds = list(map(lambda x: int(f.readline()), range(n_sp)))

        read = False
        for line in f:
            if read and not line == '\n':
                nums = line.split(' ')
                nums_clean = nums[:-1] + [nums[-1][:-2]]
                dlist.append(list(map(float, filter(bool, nums_clean))))
            if line == "\n":
                read = True

        n_trees = len(dlist)

    # make distace matrix array
    d = np.array(list(map(vector2matrix, dlist)))

    # get max distances in each tree
    max_dist = max(map(np.amax, d))

    # get max squared difference in distances between individuals
    # in lieu of epsilon

    def pairs(x):
        assert len(x) > 1
        return list(itertools.combinations(x, 2))
    def inds_by_sp(sp):
        start = sum(n_inds[:sp])
        return range(start, start + n_inds[sp])
    def dist(tup):
        i, j, n, m = itertools.chain.from_iterable(tup)
        return math.pow(d[i][n][m] - d[j][n][m], 2)

    # pairs of trees
    tree_indices = pairs(range(n_trees))
    # 'list' of lists of pairs of individuals in the same species
    nested_inds = map(compose(pairs, inds_by_sp), range(n_sp))
    # flatten the 'list' to create a 'list' of pairs
    ind_indices = itertools.chain.from_iterable(nested_inds)

    max_dist = max(map(dist, itertools.product(tree_indices, ind_indices)))


