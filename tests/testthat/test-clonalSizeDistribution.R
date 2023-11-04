# test script for clonalSizeDistribution.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalSizeDistribution works", {
  expect_doppelganger(
    "clonalSizeDistribution_plot",
    clonalSizeDistribution(combined, 
                            cloneCall = "aa", 
                            method = "ward.D2")
  )
  expect_equal(
    clonalSizeDistribution(combined, 
                          cloneCall = "aa", 
                          method = "ward.D2", 
                          exportTable = TRUE),
    getdata("visualizations", "clonalSizeDistribution_exportTable")
    )
  
})