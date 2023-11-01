# test script for clonalScatter.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalScatter works", {
  set.seed(42)
  expect_doppelganger(
    "clonalScatter_plot",
    clonalScatter(combined, 
                  cloneCall = "gene", 
                  y.axis = "P17B", 
                  x.axis = "P17L")
  )
##TODO add exportTable
})