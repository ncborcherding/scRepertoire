# test script for clonalOverlap.R - testcases are NOT comprehensive!

test_that("clonalOverlap works", {

  combined <- getCombined()

  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
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