timeit <- function(f, ...){
	ptm <- proc.time()
	x <- f(...)
	proc.time() - ptm
}