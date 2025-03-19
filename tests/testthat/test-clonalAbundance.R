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
  
  expect_doppelganger(
    "clonalAbundance_order_plot",
    clonalAbundance(getCombined(),
                    order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  combined <- getCombined()
  combined <- addVariable(combined, 
                          variable.name = "Type", 
                          variables = rep(c("B", "L"), 4))
  expect_doppelganger(
    "clonalAbundance_group_plot",
    clonalAbundance(combined, group.by = "Type")
  )
  
})