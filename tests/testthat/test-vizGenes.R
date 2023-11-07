# test script for vizGenes.R - testcases are NOT comprehensive!

test_that("vizGenes works", {

  combined <- getCombined()

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
