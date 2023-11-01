# test script for clonalProportion.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalProportion works", {
  expect_doppelganger(
    "clonalProportion_plot",
    clonalProportion(combined, 
                    cloneCall = "gene")
  )
  expect_equal(
    clonalProportion(combined, 
                    cloneCall = "gene", 
                    exportTable = TRUE),
      getdata("visualizations", "clonalProportion_exportTable")
  )
})