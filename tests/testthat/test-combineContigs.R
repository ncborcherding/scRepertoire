# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works", {
	data("contig_list")
	
	trial1 <- combineTCR(
		input.data  = lapply(contig_list[1:3], head),
		samples = c("P17B", "P17L", "P18B")
	)
	
	expected1 <- getdata("combineContigs", "combineTCR_list_expected")
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

	expect_identical(trial3, getdata("combineContigs", "combineTCR_list_filterMulti"))
	
	trial4 <- combineTCR(
	  input.data  = lapply(contig_list[1:3], head),
	  samples = c("P17B", "P17L", "P18B"),
	  removeNA = TRUE
	)
	
	expect_identical(trial4, getdata("combineContigs", "combineTCR_list_removeNA"))
}) 

# TODO combineTCR & combineBCR (need more edge cases, different args, errors, etc.)


test_that("combineBCR works", {
  
  BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
  trial1 <- combineBCR(BCR, 
                    samples = "Patient1")
  
  expect_identical(trial1, getdata("combineContigs", "combineBCR_list_expected"))
  
})