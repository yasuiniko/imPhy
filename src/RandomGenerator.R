#!/usr/bin/Rscript
'Usage: 
    RandomGenerator.R <individuals>... -e=<genes> [-p=<prob> -o=<outfile>]
    RandomGenerator.R <infile> -s=<species> [-p=<prob> -o=<outfile>]

options:
 -h --help                          Show this screen.

 -e <genes>, --genes <genes>        Number of genes.
                                    [default: 1]

 <infile>                           Filepath if using existing trees.
                                    [default: NULL]

 <individuals>                      Number of individuals per species.
                                    [default: 4]

 -o <outfile>, --outfile <outfile>  Root of filepath of outfile.
                                    (script will output outfile.txt and
                                    outfile_true.txt).
                                    [default: test]

 -p <prob>, --prob <prob>           Probability of a leaf going missing.
                                    [default: 0.2]

 -s <species>, --species <species>  Number of species.
                                    [default: 5]

Either reads a .nex file, or creates num_genes trees and randomly 
removes their leaves with probability p.' -> doc

library('docopt')
library('ape')
source('tools.R')
# Load Libraries

# Some useful functions
drop_n_tips <- function(tree, n){
    indices <- sample(as.numeric(length(tree$tip.label)), n)
    drop.tip(tree, tip=indices)
}

drop_random_tips <- function(tree, rand_fun){
    tips <- sample(as.numeric(length(tree$tip.label)), rand_fun())
    drop.tip(tree, tip=tips)
}

drop_one <- function(tree){
    dist_matrix <- cophenetic(tree)

    # select a leaf to ignore (JUST ONE FOR NOW)
    sel <- sample(as.numeric(length(tree$tip.label)), 1)

    # write the actual distance
    distances <- as.vector(as.dist(dist_matrix))

    # ignore the leaf
    dist_matrix[sel,] <- -1
    dist_matrix[,sel] <- -1

    # return both vectorized distance matrices
    list(distances, as.vector(as.dist(dist_matrix)))
}

write_matrices <- function(distances){
    actualDist <- distances[1]
    alteredDist <- distances[2]

    # write actual distance matrix
    write(paste(c(format(actualDist, digits=2)), collapse=" "), 
          fnameTRUE,
          append=TRUE)
    # write altered distance matrix
    write(paste(c(format(alteredDist, digits=2)), collapse=" "),
          fname,
          append=TRUE)
}

# Quality of Life setup
options(digits=2, error=traceback)
opts <- docopt(doc)

# Set output files -- fname=the output in the format needed, fnameTRUE=the true distances
fname <- paste(opts$outfile, '.txt', sep='')
fnameTRUE <- paste(opts$outfile, '_true.txt', sep='')
prob_missing <- assert_valid_prob(as.numeric(opts$prob))

if (!is.null(opts$infile)) {
    # get trees
    trees <- read.nexus(opts$infile)

    n_species <- as.numeric(opts$species)
    n_genes <- length(trees)
    n_leaves <- length(trees[[1]]$tip.label)
    n_individuals <- replicate(n_species, n_leaves / n_species)

    # V1 is species, V2...Vn is individuals. NA means no individual in that slot.
    si_link <- read.csv(opts$infile, header=FALSE, nrows=n_species)

} else {
    n_species <- length(opts$individuals)
    n_individuals <- as.numeric(opts$individuals)
    n_genes <- as.numeric(opts$genes)
    n_leaves <- sum(n_individuals)

    trees <- lapply(vector(,n_genes), function(x) rcoal(n_leaves))
}

# Write the number of species, number of leaves to the file
invisible(write(c(n_species, n_leaves), fname, append=TRUE))

# write the number of individuals to file
invisible(lapply(n_individuals, write, file=fname, append=TRUE))
invisible(write("\n", fname, append=TRUE))

invisible(lapply(trees, function(x) write_matrices(drop_one(x))))
