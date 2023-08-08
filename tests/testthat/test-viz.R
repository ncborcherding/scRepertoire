# test script for viz.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(
	contig_list, 
	samples = c("PY", "PY", "PX", "PX", "PZ","PZ"), 
	ID = c("P", "T", "P", "T", "P", "T")
)

test_that("quantContig works", {
	expect_doppelganger(
		"quantContig_scaled_plot", quantContig(combined, scale = TRUE)
	)
	
	# TODO test for only a single sample
	# TODO test the exportTable param
})

test_that("abundanceContig works", {
	expect_doppelganger(
		"abundanceContig_scaled_plot", abundanceContig(combined, scale = FALSE)
	)
})

test_that("lengthContig works", {
	expect_doppelganger(
		"lengthContig_both_chain_plot", lengthContig(combined, chain = "both") 
	)
})

test_that("compareClonotypes works", {
	expect_doppelganger(
		"compareClonotypes_vignette_plot",
		compareClonotypes(
			combined, 
			numbers = 10, 
			samples = c("PX_P", "PX_T"), 
			cloneCall="aa", 
			graph = "alluvial"
		)
	)
})

test_that("vizGenes works", {
	expect_doppelganger(
		"vizGenes_bar_vignette_plot",
		vizGenes(
			combined,
			gene = "V", 
			chain = "TRB", 
			plot = "bar", 
			order = "variance", 
			scale = TRUE
		)
	)
	
	expect_doppelganger(
		"vizGenes_heatmap_vignette_plot",
		vizGenes(
			combined[c(1,3,5)], 
			gene = "V", 
			chain = "TRB", 
			y.axis = "J", 
			plot = "heatmap", 
			scale = TRUE, 
			order = "gene"
		)
	)
})

# TODO scatterClonotype
# TODO clonesizeDistribution
# TODO makingLodes
# TODO vizGenes
