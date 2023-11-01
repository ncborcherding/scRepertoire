# test script for percentAA.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("percentAA works", {
  set.seed(42)
  expect_doppelganger(
    "percentAA_plot",
    percentAA(combined, 
              chain = "TRB", 
              aa.length = 20)
  )
  set.seed(42)
  expect_equal(
    percentAA(combined, 
              chain = "TRB", 
              aa.length = 20,
              exportTable = TRUE),
      getdata("visualizations", "percentAA_exportTable")
  )
})