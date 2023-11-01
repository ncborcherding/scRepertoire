# test script for clonalDiversity.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalDiversity works", {
  set.seed(42)
  expect_doppelganger(
    "clonalDiversity_plot",
    clonalDiversity(combined, 
                    cloneCall = "gene", 
                    n.boots = 1)
  )
  set.seed(42)
  expect_equal(
      clonalDiversity(combined, 
                      cloneCall = "gene", 
                      skip.boots = TRUE, 
                      metrics = c("norm.entropy", "shannon"),
                      exportTable = TRUE),
      getdata("visualizations", "clonalDiversity_exportTable")
  )
})