library('ape')
library('phytools')
# Load Libraries

# Set output files -- fname=the output in the format needed, fnameTRUE=the true distances
fname <- "test.txt"
fnameTRUE <- "trueTest.txt"

# Parameters: numSpecies, numLeafs, probability of missing leaf
numSpecies <- 5
numLeafs <- 4
prob_missing <- 0.2

# Output set to 2 decimal digits
options(digits=2)

# Write the number of species, number of leafs in the file
write(c(numSpecies,numLeafs), fname, append=TRUE)

# Now I have the population of every species to be up to 10: THIS CAN BE CHANGED
poss_Nums <- 1:10

# Creating Random numbers of individuals for each species
numIndividuals <- c();
for (j in 1:numSpecies){
    numIndividuals <- c(numIndividuals, sample.int(10, size = 1, replace = TRUE));
}

# Writing the numbers of individuals in file
for (j in 1:numSpecies){
    write(numIndividuals[j], fname, append=TRUE)
}
write("\n", fname, append=TRUE)


# Creating a Sample to choose from for leafs missing
sam <- 1:numLeafs

# For every species
for (j in 1:numSpecies){
    # For every individual in that species
    for (i in 1:numIndividuals[j]){
        # Create random tree with numLeafs leafs
        tree <- rcoal(numLeafs)
        # Calculate the distances in vectorized form
        dist <- cophenetic(tree)
        # Randomly select some trees to be missing a leaf
        random <- runif(1, 0, 1)
        if (random < prob_missing){
            # Randomly select which leaf is missing: NOW ONLY ONE LEAF MISSING
            sel <- sample(sam, 1)
        }
        else{
            sel <- -1
        }
        # Create vectors for actual and altered Distances
        actualDist <- c()
        alteredDist <- c()
        for (k in 1:(numLeafs-1)){
            for (l in (k+1):numLeafs){
                x1 <- paste("t",k,sep="")
                y1 <- paste("t",l,sep="")
                # Write the actual distances in the format needed
                actualDist <- c(actualDist, dist[x1,y1])
                # If the leaf is missing, set its altered distance to -1
                if (sel==k || sel==l){
                    alteredDist <- c(alteredDist, -1)
                }
                else{
                    alteredDist <- c(alteredDist, dist[x1,y1])
                }
            }
        }
        # Write altered and actual distances in the format of the file
    write(paste(c(format(actualDist, digits=2)), collapse=" "), fnameTRUE ,append=TRUE)
    
    write(paste(c(format(alteredDist, digits=2)), collapse=" "), fname ,append=TRUE)
    }
    # Leave an empty line between species (for easy visual separation)
    write("\n", fnameTRUE, append=TRUE)
    write("\n", fname, append=TRUE)
   
}

