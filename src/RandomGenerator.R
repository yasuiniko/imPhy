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

# Some useful functions
drop_tips <- function(tree, n_tips_to_drop){

    # organize distance matrix
    d <- cophenetic(tree)
    dist_matrix <- d[order(rownames(d)), order(colnames(d))]

    # select a leaf to ignore
    # for random: n_tips_to_drop <- rbinom(1, num_tips, prob_missing)
    sel <- sample(n_leaves, n_tips_to_drop(n_leaves))

    # save the actual distances
    distances <- as.dist(dist_matrix)

    # ignore the leaf
    dist_matrix[sel,] <- -1
    dist_matrix[,sel] <- -1

    # return both vectorized distance matrices
    list(distances, as.dist(dist_matrix))
}

write_matrices <- function(distances){
    actualDist <- distances[1]
    alteredDist <- distances[2]

    # write actual distance matrix
    write(paste(sapply(actualDist, two_digits), collapse=' '),
          fnameTRUE,
          append=TRUE)
    # write altered distance matrix
    write(paste(sapply(alteredDist, two_digits), collapse=' '),
          fname,
          append=TRUE)
}

# Quality of Life setup
options(digits=2, error=traceback)
opts <- docopt(doc)

# Set output files -- fname=the output in the format needed, fnameTRUE=the true distances
fname <- paste0(opts$outfile, '.txt')
fnameTRUE <- paste0(opts$outfile, '_true.txt')
prob_missing <- assert_valid_prob(as.numeric(opts$prob))
n_tips_to_drop <- function(max_tips) rbinom(1, max_tips, prob_missing)

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

invisible(lapply(trees, function(x) write_matrices(drop_tips(x, n_tips_to_drop))))
