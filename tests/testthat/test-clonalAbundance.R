# test script for clonalAbundance.R - testcases are NOT comprehensive!

test_that("clonalAbundance works", {
  expect_doppelganger(
    "clonalAbundance_scaled_plot",
    clonalAbundance(getCombined(), scale = FALSE)
  )
})
