# test script for positionalProperty.R - testcases are NOT comprehensive!

test_that("positionalProperty works", {

  set.seed(42)
  expect_doppelganger(
    "positionalEntropy_TRB_plot",
     positionalProperty(getCombined(), 
                        chain = "TRB", 
                        aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_TRB_order_plot",
    positionalProperty(getCombined(), 
                       chain = "TRB", 
                       aa.length = 20,
                       order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  expect_doppelganger(
    "positionalEntropy_TRA_plot",
    positionalProperty(getCombined(), 
                       chain = "TRA", 
                       aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_Kidera_plot",
    positionalProperty(getCombined(), 
                       chain = "TRB", 
                       method = "Kidera",
                       aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_stScales_plot",
    positionalProperty(getCombined(), 
                       chain = "TRB", 
                       method = "stScales",
                       aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_tScales_plot",
    positionalProperty(getCombined(), 
                       chain = "TRB", 
                       method = "tScales",
                       aa.length = 20)
  )
  
  expect_doppelganger(
    "positionalEntropy_VHSE_plot",
    positionalProperty(getCombined(), 
                       chain = "TRB", 
                       method = "VHSE",
                       aa.length = 20)
  )
})