library('ape')

"
USAGE: To drop 4 random tips from each phylo (tree) in the multiPhylo

	my_multiPhylo <- read.nexus('treefile.nex')

	source('path/to/remove_random_tips.R')
	dropped_multiPhylo <- drop_tips(my_multiPhylo, 4)
"


.drop_indices <- function(tree, num_to_drop){
	sample(as.numeric(length(tree$tip.label)), num_to_drop)
}

.drop_tip <- function(tree, n){
	drop.tip(tree, tip=.drop_indices(tree, num_to_drop=n))
}

drop_tips <- function(trees, num_to_drop){
	lapply(trees, .drop_tip, n=num_to_drop)
}