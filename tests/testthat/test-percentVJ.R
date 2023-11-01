# test script for percentVJ.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("percentVJ works", {
  set.seed(42)
  expect_doppelganger(
    "percentVJ_plot",
    percentVJ(combined, 
              chain = "TRB")
  )
  set.seed(42)
  expect_equal(
    percentVJ(combined, 
              chain = "TRB", 
              exportTable = TRUE),
      getdata("visualizations", "percentVJ_exportTable")
  )
})