# test script for clonalHomeostasis.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalDiversity works", {
  expect_doppelganger(
    "clonalHomeostasis_plot",
    clonalHomeostasis(combined, 
                    cloneCall = "gene")
  )
  expect_equal(
    clonalHomeostasis(combined, 
                      cloneCall = "gene", 
                      exportTable = TRUE),
      getdata("visualizations", "clonalHomeostasis_exportTable")
  )
})