# test script for clonalQuant.R - testcases are NOT comprehensive!

test_that("clonalQuant works", {

	combined <- getCombined()
	combined <- addVariable(combined, 
	                        variable.name = "Type", 
	                        variables = rep(c("B", "L"), 4))

	expect_doppelganger(
		"clonalQuant_scaled_plot", 
		clonalQuant(combined, 
		            scale = TRUE)
	)
	
	expect_doppelganger(
	  "clonalQuant_order_plot", 
	  clonalQuant(combined, 
	              scale = TRUE,
	              order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
	)
	
	expect_doppelganger(
	  "clonalQuant_unscaled_plot",
  	clonalQuant(combined, 
  	            group.by = "Type",
  	            scale = FALSE)
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
