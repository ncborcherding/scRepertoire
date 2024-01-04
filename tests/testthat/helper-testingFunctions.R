getdata <- function(dir, name) {
	readRDS(paste("testdata/", dir, "/", name, ".rds", sep = "")) # could move testdata 1 dir lvl up nstead
}

# shortcut get the combined contigs that are used in most testcases
getCombined <- function() {
	getdata("combineContigs", "combined")
}
