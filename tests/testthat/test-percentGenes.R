# test script for percentGenes.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("percentGenes works", {
  set.seed(42)
  expect_doppelganger(
    "percentGenes_plot",
    percentGenes(combined, 
                 chain = "TRB", 
                 gene = "V")
  )
  set.seed(42)
  expect_equal(
    percentGenes(combined, 
              chain = "TRB", 
              gene = "V",
              exportTable = TRUE),
      getdata("visualizations", "percentGenes_exportTable")
  )
})