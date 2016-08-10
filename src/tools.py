from functools import reduce
import itertools
import multiprocessing
import time

def timeit(f, s):
    bigs = s[0].upper() + s[1:]
    smalls = s[0].lower() + s[1:]
    print("{}...".format(bigs))
    t = time.time()
    x = f()
    print("Done {} in {}s.\n".format(smalls, time.time() - t))
    return x

def flatten(*lsts):
    return itertools.chain(*lsts)

def iterflatten(iterable):
    return itertools.chain.from_iterable(iterable)

def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    """
    Taken from http://stackoverflow.com/revisions/16071616/9
    """
    def fun(f, q_in, q_out):
        while True:
            i, x = q_in.get()
            if i is None:
                break
            q_out.put((i, f(x)))

    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
            for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]

gene =           'd{}_g{}_i{}_n{}_s{}_e{}.nex'
batch_analyze =  "d{0}_g{2}_i{3}_n{5}_s{7}"
batch_general =  "d{}_g{}_i{}_n{}_s{}"
fileroot =       "d{0}_g{2}_i{3}_n{5}_s{7}_e{1}_m{4}_p{6}"

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

        # species_taxa = [taxa[i:i+k] for i in range(0, len(taxa), k)]

        return all(map(has_species_subtree, species_taxa))

    num_separated = sum(map(is_separated_tree, treelist))
    total = len(treelist)

    return num_separated / total