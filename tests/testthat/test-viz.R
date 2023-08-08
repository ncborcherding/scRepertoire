# test script for viz.R - testcases are NOT comprehensive!
# TODO all functions need to be tested for a single sample - however, many don't
# really work properly for single samples

# testdata: (assumes combineTCR works)
combined <- combineTCR(
	contig_list, 
	samples = c("PY", "PY", "PX", "PX", "PZ","PZ"), 
	ID = c("P", "T", "P", "T", "P", "T")
)

single_contig <- combineTCR(contig_list[[1]])

single_contig_with__sample <- combineTCR(
	contig_list[[1]], samples = "PX", ID = "P"
)

test_that("quantContig works", {
	expect_doppelganger(
		"quantContig_scaled_plot", quantContig(combined, scale = TRUE)
	)
	expect_equal(
		quantContig(combined, scale = TRUE, exportTable = TRUE),
		data.frame(
			"contigs" = c(2712L, 1585L, 823L, 918L, 1143L, 768L),
			"values" = c("PY_P", "PY_T", "PX_P", "PX_T", "PZ_P", "PZ_T"), 
			"total" = c(3208L, 3119L, 1068L, 1678L, 1434L, 2768L),
			"scaled" = c(
				84.5386533665835, 50.8175697338891, 77.059925093633,
				54.7079856972586, 79.7071129707113, 27.7456647398844
			)
		)
	)
	
	expect_doppelganger(
		"quantContig_single_sample_plot", quantContig(single_contig)
	)
	expect_identical(
		quantContig(single_contig, exportTable = TRUE),
		data.frame("contigs" = 2712L, "total" = 3208L)
	)
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
		"compareClonotypes_alluvial_plot",
		compareClonotypes(
			combined, 
			numbers = 10, 
			samples = c("PX_P", "PX_T"), 
			cloneCall="aa", 
			graph = "alluvial"
		)
	)
	
	expect_doppelganger(
		"compareClonotypes_area_plot",
		compareClonotypes(
			combined, 
			numbers = 10, 
			samples = c("PX_P", "PX_T"), 
			cloneCall="aa", 
			graph = "area"
		)
	)
})

test_that("scatterClonotype works", {
	expect_doppelganger(
		"scatterClonotype_vignette_plot",
		scatterClonotype(
			combined, 
			cloneCall ="gene", 
			x.axis = "PY_P", 
			y.axis = "PY_T",
			dot.size = "total",
			graph = "proportion",
			seed = 42
		)
	)
	
	# TODO test the exportTable arg
})

test_that("clonesizeDistribution works", {
	expect_doppelganger(
		"clonesizeDistribution_vignette_plot",
		clonesizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
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
