# test script for clonalLength.R - testcases are NOT comprehensive!

test_that("clonalLength works", {
  expect_doppelganger(
    "clonalLength_both_chain_plot", clonalLength(getCombined(), chain = "both") 
  )
})

