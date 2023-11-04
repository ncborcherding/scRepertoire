# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works", {
	data("contig_list")
	
	trial1 <- combineTCR(
		input.data  = lapply(contig_list[1:3], head),
		samples = c("P17B", "P17L", "P18B")
	)
	
	expected1 <- readRDS("testdata/combineContigs/combineTCR_list_expected.rds")
	
	expect_identical(trial1, expected1)
	
	trial2 <- combineTCR(
		head(contig_list[[1]]), samples = "P17B"
	)[[1]]
	expected2 <- expected1[[1]]
	
	expect_identical(trial2, expected2)
	
	trial3 <- combineTCR(
	  input.data  = lapply(contig_list[1:3], head),
	  samples = c("P17B", "P17L", "P18B"),
	  filterMulti = TRUE
	)
	expected3 <- readRDS("testdata/combineContigs/combineTCR_list_filterMulti.rds")
	
	expect_identical(trial3, expected3)
	
	trial4 <- combineTCR(
	  input.data  = lapply(contig_list[1:3], head),
	  samples = c("P17B", "P17L", "P18B"),
	  removeNA = TRUE
	)
	expected4 <- readRDS("testdata/combineContigs/combineTCR_list_removeNA.rds")
	
	expect_identical(trial4, expected4)
	
}) 

# TODO combineTCR (need more edge cases, different args, errors, etc.)
# TODO combineBCR
