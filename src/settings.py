from itertools import product

class Settings:

    """
    Testing variables.
    """
    t_c = [1, 10, 20]                 # c ratio
    t_genes = [3,4,5]                 # number of genes
    t_inds = [8]                      # number of individuals per species
    t_methods = [1]                   # imputation methods to use
    t_probs = [4, 8, 16]              # leaf dropping probabilities/denominators
    t_species = [2, 4, 6]             # number of species 
    t_trees = [3]                     # number of species trees
    t_pop_sizes = [10000]             # effective population size
    t_dists = ['norm']                # distances to use

    t_c = list(map(float, t_c))
    t_depths = list(set(map(lambda x: int(x[0]*x[1]), product(t_c,
                                                              t_pop_sizes))))
    t_flow_dict = {"generate":False,  # generate trees
                   "drop":True,      # drop leaves
                   "impute":False,    # impute missing leaves
                   "analyze":False}    # analyze batches

    """
    Experiment variables.
    """
    e_c = [1, 4, 8, 12, 16, 20]       # c ratio
    e_genes = [10, 20, 30, 200, 500]  # number of genes
    e_inds = [8]                      # number of individuals per species
    e_methods = [1]                   # imputation methods to use
    e_probs = [4, 8, 16]              # leaf dropping probabilities/denominators
    e_species = [2, 4, 6, 8]          # number of species 
    e_trees = [3]                     # number of species trees
    e_pop_sizes = [10000]              # effective population size
    e_dists = ['norm']   # distances to use

    # convert c list from integers to floats
    e_c = list(map(float, e_c))
    # depth (in generations) is equal to c*pop_size, for each 
    # combination of values of c and pop_size. If you prefer, 
    # replace the following line with a hard-coded depth list, as 
    # seen above in the definition for c. 
    e_depths = list(set(map(lambda x: int(x[0]*x[1]), product(e_c, 
                                                              e_pop_sizes))))

    # Options to set 
    e_flow_dict = {"generate":False,  # generate trees
                 "drop":False,      # drop leaves
                 "impute":False,    # impute missing leaves
                 "analyze":True}    # analyze batches

    def __init__(self,
                 c,
                 genes,
                 inds,
                 methods,
                 probs,
                 species,
                 trees,
                 pop_sizes,
                 dists,
                 flow_dict,
                 depths=[]):

        self.c = list(map(float, c))
        self.genes = genes
        self.inds = inds
        self.methods = methods
        self.probs = probs
        self.species = species
        self.trees = trees
        self.pop_sizes = pop_sizes
        self.dists = dists

        if depths:
            self.depths = depths
        else:
            self.depths = list(set(map(lambda x: int(x[0]*x[1]), 
                                       product(self.c, self.pop_sizes))))

        self.flow_dict = flow_dict

    @staticmethod
    def get_e():
        return Settings(Settings.e_c,
                        Settings.e_genes,
                        Settings.e_inds,
                        Settings.e_methods,
                        Settings.e_probs,
                        Settings.e_species,
                        Settings.e_trees,
                        Settings.e_pop_sizes,
                        Settings.e_dists,
                        Settings.e_flow_dict)

    @staticmethod
    def get_t():
        return Settings(Settings.t_c,
                        Settings.t_genes,
                        Settings.t_inds,
                        Settings.t_methods,
                        Settings.t_probs,
                        Settings.t_species,
                        Settings.t_trees,
                        Settings.t_pop_sizes,
                        Settings.t_dists,
                        Settings.t_flow_dict)

