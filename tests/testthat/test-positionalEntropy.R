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
    "positionalEntropy_TRB_order_plot",
    positionalEntropy(getCombined(), 
                      chain = "TRB", 
                      aa.length = 20,
                      order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  expect_doppelganger(
    "positionalEntropy_TRA_plot",
    positionalEntropy(getCombined(), 
                      chain = "TRA", 
                      aa.length = 20)
  )
})