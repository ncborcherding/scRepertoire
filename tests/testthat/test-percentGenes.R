# test script for percentGenes.R - testcases are NOT comprehensive!

test_that("percentGenes works", {

  combined <- getCombined()

  set.seed(42)
  expect_doppelganger(
    "percentGenes_plot",
    percentGenes(combined, 
                 chain = "TRB", 
                 gene = "V")
  )
  set.seed(42)
  expect_equal(
    percentGenes(combined, 
              chain = "TRB", 
              gene = "V",
              exportTable = TRUE),
      getdata("visualizations", "percentGenes_exportTable")
  )
})