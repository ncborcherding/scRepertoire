# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works", {
	data("contig_list")
	trial <- combineTCR(
		df = lapply(contig_list[1:3], head),
		samples = c("PY", "PY", "PX"),
		ID = c("P", "T", "P")
	)
	expect_identical(
		trial,
		readRDS("testdata/combineContigs/combineTCR_list_expected.rds")
	)
}) 

# TODO combineTCR (need more edge cases, different args, errors, etc.)
# TODO combineBCR
# TODO lvCompare
