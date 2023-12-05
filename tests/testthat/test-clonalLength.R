# test script for clonalLength.R - testcases are NOT comprehensive!

test_that("clonalLength works", {
  expect_doppelganger(
    "clonalLength_both_chain_plot", 
    clonalLength(getCombined(), chain = "both") 
  )
  
  expect_doppelganger(
    "clonalLength_scaled_plot", 
    clonalLength(getCombined(), 
                 scale = TRUE) 
  )
  
  combined <- addVariable(getCombined(), 
                          variable.name = "Type", 
                          variables = rep(c("B", "L"), 4))
  expect_doppelganger(
    "clonalLength_groupby_plot",
    clonalLength(combined, 
                    cloneCall = "nt",
                    group.by = "Type")
  )
})

