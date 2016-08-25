"""
This file is part of imPhy, a pipeline for evaluating the quality of
phylogenetic imputation software.
Copyright © 2016 Niko Yasui, Chrysafis Vogiatzis

imPhy uses GTP, which is Copyright © 2008, 2009  Megan Owen, Scott Provan

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import decimal
from functools import reduce
import gzip
import itertools
import logging
import multiprocessing
import os
import shutil
import subprocess
import sys
import threading
import time


########### CONSTANTS ###########
gene =           'd{}_g{}_i{}_n{}_s{}_e{}.nex'
batch_analyze =  "d{0}_g{2}_i{3}_n{5}_s{7}"
batch_general =  "d{}_g{}_i{}_n{}_s{}"
fileroot =       "d{0}_g{2}_i{3}_n{5}_s{7}_e{1}_m{4}_p{6}"
ext = ".txt.gz"
ext_len = len(ext)

########### FUNCTIONS ###########

def gunzip_to(dest_root, src):
    dest = os.path.join(dest_root, os.path.basename(src)[:-3])
    with gzip.open(src, 'rt') as f_in:
        with open(dest, 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)

def gzip_to(dest_root, src):
    dest = os.path.join(dest_root, os.path.basename(src)+'.gz')
    with open(src, 'rt') as f_in:
        with gzip.open(dest, 'wt') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(src)

def to_decimal(x):
    return decimal.Decimal(str(x))

def timeit(f, s, logger=logging.getLogger()):
    big_s = s[0].upper() + s[1:]
    small_s = s[0].lower() + s[1:]

    logger.info("{}...".format(big_s))

    t = time.time()
    x = f()

    logger.info("Done {} in {}s.".format(small_s, time.time() - t))

    return x

def flatten(*lsts):
    return itertools.chain(*lsts)

def iterflatten(iterable):
    return itertools.chain.from_iterable(iterable)

def get_output(cargs, logger=logging.getLogger(), ignore_error=False):
    """
    Calls subprocess.check_output on the args and logs the output.
    """

    try:
        output = str(subprocess.check_output(cargs), "utf-8")
        if output:
            logger.debug("Output of {}: {}".format(cargs, output))

    except subprocess.CalledProcessError as e:
        if ignore_error:
            output = "error"
        logger.warning("CalledProcessError when calling {}.".format(cargs),
                       exc_info=True)

    return output

class MultiLogger:
    def __init__(self, logpath):
        # self.q_log = multiprocessing.Queue()

        root = logging.getLogger()
        root.setLevel(logging.DEBUG)
        
        logfile_h = logging.FileHandler(logpath, 'w')
        log_format = '%(asctime)s %(processName)-10s'
        log_format += ' %(name)s %(levelname)-8s %(message)s'
        logfile_h.setFormatter(logging.Formatter(log_format))
        logfile_h.setLevel(logging.DEBUG)
        root.addHandler(logfile_h)

        std_out_h = logging.StreamHandler(sys.stdout)
        std_out_h.setLevel(logging.INFO)
        root.addHandler(std_out_h)

    def getLogger(self, name=""):
        return logging.getLogger(name)

def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    """
    Code adapted from http://stackoverflow.com/revisions/16071616/9
    """

    def worker(f, q_in, q_out):
        child_logger = logging.getLogger(__name__)

        while True:
            i, x = q_in.get()
            if i is None:
                break
            q_out.put((i, f(x)))


    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    worker_args = (f, q_in, q_out)
    new_worker = lambda: multiprocessing.Process(target=worker, 
                                                 args=worker_args)

    proc = [new_worker() for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]

# def parmap(f, X, nprocs=multiprocessing.cpu_count()):
#     """
#     Taken from http://stackoverflow.com/revisions/16071616/9
#     Here as a backup in case my logging doesn't work.
#     """
#     def fun(f, q_in, q_out):
#         while True:
#             i, x = q_in.get()
#             if i is None:
#                 break
#             q_out.put((i, f(x)))

#     q_in = multiprocessing.Queue(1)
#     q_out = multiprocessing.Queue()

#     proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
#             for _ in range(nprocs)]
#     for p in proc:
#         p.daemon = True
#         p.start()

#     sent = [q_in.put((i, x)) for i, x in enumerate(X)]
#     [q_in.put((None, None)) for _ in range(nprocs)]
#     res = [q_out.get() for _ in range(len(sent))]

#     [p.join() for p in proc]

#     return [x for i, x in sorted(res)]

def prop_separated_trees(treelist, get_species=lambda taxon: taxon.label[0]):
    """
    Calculates the proportion of trees for which each tree has a subtree
    containing all the members of a species and none of the members of 
    another species.
    """

    def is_separated_tree(tree):
        """
        Determines if each species in the tree can be contained in a 
        subtree which contains only that species.
        """

        def has_species_subtree(sp_taxa):
            """
            Check if the subtree whose leaves contain all taxa in 
            sp_taxa has leaves with taxa outside of sp_taxa.
            """

            # initialize labels and node
            labels = [n.label for n in sp_taxa]
            node = tree.find_node_with_taxon_label(labels[0])
            nlabs = [n.label for n in node.leaf_nodes()]

            # loop until all labels are in node's subtree
            while not all(lab in nlabs for lab in labels):
                node = node.parent_node
                nlabs = [n.taxon.label for n in node.leaf_nodes()]

            # check if there are any outspecies leaves in the subtree
            return len(nlabs) == len(labels)

        # find the number of individuals per species
        n_sp = len(set(map(get_species, tree.taxon_namespace)))
        k = int(len(tree.taxon_namespace)/n_sp)

        # split taxa into groups by species
        taxa = list(sorted(tree.taxon_namespace, key=get_species))
        species_taxa = []
        while taxa:
            curr_species = get_species(taxa[0])
            in_species = lambda t: get_species(t) == curr_species
            species = list(filter(in_species, taxa))
            species_taxa.append(species)
            list(map(taxa.remove, species))

        is_separated = all(map(has_species_subtree, species_taxa))

        return is_separated

    num_separated = sum(map(is_separated_tree, treelist))
    total = len(treelist)

    return num_separated / total