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
  
  combined <- getCombined()
  combined <- addVariable(getCombined(), 
                          variable.name = "Type", 
                          variables = rep(c("B", "L"), 4))
  expect_doppelganger(
    "clonalAbundance_group_plot",
    clonalAbundance(combined, group.by = "Type")
  )
  
})
#TODO Add grouping plot