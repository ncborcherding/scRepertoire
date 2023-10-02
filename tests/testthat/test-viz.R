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
		clonalQuant(combined, scale = TRUE, exportTable = TRUE),
		data.frame(
			"contigs" = c(745L, 2117L, 1254L, 1202L, 5544L, 1619L, 6087L,  192L),
			"values" = c("P17B", "P17L", "P18B", "P18L", "P19B","P19L", "P20B", "P20L"), 
			"total" = c(2805L, 2893L, 1328L, 1278L, 6942L, 2747L, 8991L,  201L),
			"scaled" = c(
			  26.55971, 73.17663, 94.42771, 94.05321, 79.86171, 58.93702, 67.70103, 95.52239
			)
		)
	)
	
	expect_doppelganger(
		"clonalQuant_single_sample_plot", clonalQuant(single_contig)
	)
	expect_identical(
		clonalQuant(single_contig, exportTable = TRUE),
		data.frame("contigs" = 745L, "total" = 2805L)
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
			numbers = 10, 
			samples = c("P17B", "P17L"), 
			cloneCall="aa", 
			graph = "alluvial"
		)
	)
	
	expect_doppelganger(
		"clonalCompare_area_plot",
		clonalCompare(
			combined, 
			numbers = 10, 
			samples = c("P17B", "P17L"), 
			cloneCall="aa", 
			graph = "area"
		)
	)
})

test_that("clonalScatter works", {
	expect_doppelganger(
		"clonalScatter_vignette_plot",
		clonalScatter(
			combined, 
			cloneCall ="gene", 
			x.axis = "P17B", 
			y.axis = "P17L",
			dot.size = "total",
			graph = "proportion",
			seed = 42
		)
	)
	
	# TODO test the exportTable arg
})

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
