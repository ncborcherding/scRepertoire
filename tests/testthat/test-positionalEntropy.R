# test script for positionalEntropy.R - testcases are NOT comprehensive!

test_that("positionalEntropy works", {

  set.seed(42)
  expect_doppelganger(
    "positionalEntropy_TRB_plot",
     positionalEntropy(getCombined(), 
                      chain = "TRB", 
                      aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_TRA_plot",
    positionalEntropy(getCombined(), 
                      chain = "TRA", 
                      aa.length = 20)
  )
})