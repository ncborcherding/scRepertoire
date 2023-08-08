# test script for viz.R - testcases are NOT comprehensive!

test_that("quantContig works, assuming combineTCR works", {
	expect_doppelganger(
		"quantContig_mini_scaled_plot",
		quantContig(
			combineTCR(
				contig_list,
				samples = c("PY", "PY", "PX", "PX", "PZ","PZ"), 
				ID = c("P", "T", "P", "T", "P", "T")
			),
			scale = TRUE
		)
	)
	
	# TODO test for only a single sample
	# TODO test the exportTable param
})

# TODO abundanceContig
# TODO lengthContig
# TODO compareClonotypes
# TODO scatterClonotype
# TODO clonesizeDistribution
# TODO makingLodes
# TODO vizGenes
