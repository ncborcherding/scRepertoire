# test script for clonalOverlap.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalOverlap works", {
  expect_doppelganger(
    "clonalOverlap_raw_plot",
    clonalOverlap(
      combined, 
      method = "raw")
  )
  expect_equal(
    clonalOverlap(
      combined[c(1,2)], 
      method = "raw", 
      exportTable = TRUE),
    data.frame(row.names = c("P17B", "P17L"),
      "P17B" = rep(NA, 2),
      "P17L" = c(85L, NA)
    )
  )
})