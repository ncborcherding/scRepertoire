# test script for viz.R - testcases are NOT comprehensive!
# TODO all functions need to be tested for a single sample - however, many don't
# really work properly for single samples

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))

single_contig <- combineTCR(contig_list[[1]])

single_contig_with_sample <- combineTCR(
	contig_list[[1]], samples = "P17B",
)
# TODO test more cases with single_contig

test_that("clonalQuant works", {
	expect_doppelganger(
		"clonalQuant_scaled_plot", clonalQuant(combined, scale = TRUE)
	)
	expect_equal(
		clonalQuant(combined, 
		            scale = TRUE, 
		            exportTable = TRUE),
		data.frame(
			"contigs" = c(745L, 2117L, 1254L, 1202L, 5544L, 1619L, 6087L,  192L),
			"values" = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"), 
			"total" = c(2805L, 2893L, 1328L, 1278L, 6942L, 2747L, 8991L,  201L),
			"scaled" = c(
			  26.5597148, 73.1766333, 94.4277108, 94.0532081, 79.8617113, 58.9370222, 67.7010344, 95.5223881
			)
		)
	)
})

test_that("clonalAbundance works", {
	expect_doppelganger(
		"clonalAbundance_scaled_plot", clonalAbundance(combined, scale = FALSE)
	)
})

test_that("clonalLength works", {
	expect_doppelganger(
		"clonalLength_both_chain_plot", clonalLength(combined, chain = "both") 
	)
})

test_that("clonalCompare works", {
	expect_doppelganger(
		"clonalCompare_alluvial_plot",
		clonalCompare(
			combined, 
			top.clones = 10, 
			samples = c("P17B", "P17L"), 
			cloneCall="aa", 
			graph = "alluvial"
		)
	)
	
	expect_doppelganger(
		"clonalCompare_area_plot",
		clonalCompare(
		  input.data = combined, 
			top.clones  = 10, 
			samples = c("P17B", "P17L"), 
			cloneCall="aa", 
			graph = "area"
		)
	)
})

# TODO clonalScatter

# something in `clonesizeDistribution` prints "NULL" to the terminal
test_that("clonalSizeDistribution works", {
	expect_doppelganger(
		"clonalSizeDistribution_vignette_plot",
		clonalSizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
	)
})

# TODO makingLodes


test_that("vizGenes works", {
	expect_doppelganger(
		"vizGenes_bar_vignette_plot",
		vizGenes(
			combined,
			x.axis = "TRBV", 
			y.axis = NULL,
			plot = "barplot", 
			order = "variance", 
			scale = TRUE
		)
	)
	
	expect_doppelganger(
		"vizGenes_heatmap_vignette_plot",
		vizGenes(
			combined[c(1,3,5)], 
			x.axis = "TRBV",
			y.axis = "TRBJ",
			plot = "heatmap", 
			scale = TRUE, 
			order = "gene"
		)
	)
})
