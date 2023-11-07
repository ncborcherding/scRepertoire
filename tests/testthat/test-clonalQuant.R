# test script for clonalQuant.R - testcases are NOT comprehensive!

test_that("clonalQuant works", {

	combined <- getCombined()

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
