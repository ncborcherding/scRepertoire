# test script for clonalSizeDistribution.R - testcases are NOT comprehensive!

test_that("clonalSizeDistribution works", {

  combined <- getCombined()
  set.seed(42)
  expect_doppelganger(
    "clonalSizeDistribution_plot",
    clonalSizeDistribution(combined, 
                            cloneCall = "aa", 
                            method = "ward.D2")
  )
  expect_equal(
    clonalSizeDistribution(combined, 
                          cloneCall = "aa", 
                          method = "ward.D2", 
                          exportTable = TRUE),
    getdata("visualizations", "clonalSizeDistribution_exportTable"),
    tolerance = 1e-4
  )
  
})