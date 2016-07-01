
library('ape')
 
################## Parse Args ##################
# should be using a library, but haven't been
# able to get one to work

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0){
	print("Using default filename 'out.tree'.")
	outfile <- "out.tree"
} else{
	outfile <- args[1]
}
################ End Parse Args ################

# overwrite outfile if it exists
if (file.exists(outfile)) {
	invisible(file.remove(outfile))
}

# only make the function once instead of using an anonymous function
make_tree <- function(x) rcoal(10)

# make 1000 trees
y <- lapply(vector(, 1000), make_tree)

# write the trees
invisible(lapply(y, write.tree, file=outfile, append=TRUE))