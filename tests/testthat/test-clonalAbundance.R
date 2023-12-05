# test script for clonalAbundance.R - testcases are NOT comprehensive!

test_that("clonalAbundance works", {
  expect_doppelganger(
    "clonalAbundance_unscaled_plot",
    clonalAbundance(getCombined(), scale = FALSE)
  )
  
  expect_doppelganger(
    "clonalAbundance_scaled_plot",
    clonalAbundance(getCombined(), scale = TRUE)
  )
})
#TODO Add grouping plot