# test script for clonalScatter.R - testcases are NOT comprehensive!

test_that("clonalScatter works", {
  set.seed(42)
  expect_doppelganger(
    "clonalScatter_plot",
    clonalScatter(getCombined(), 
                  cloneCall = "gene", 
                  y.axis = "P17B", 
                  x.axis = "P17L")
  )
##TODO add exportTable
})
