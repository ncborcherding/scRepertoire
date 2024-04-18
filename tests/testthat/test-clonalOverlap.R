# test script for clonalOverlap.R - testcases are NOT comprehensive!

test_that("clonalOverlap works", {

  combined <- getCombined()

  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_raw_plot",
    clonalOverlap(
      combined, 
      method = "raw")
  )
  
  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_reorder_plot",
    clonalOverlap(
      combined[c(3,4,1,2,7,8,5,6)], 
      method = "raw")
  )
  
  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_morisita_plot",
    clonalOverlap(
      combined, 
      method = "morisita")
  )
  
  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_jaccard_plot",
    clonalOverlap(
      combined, 
      chain = "TRA",
      method = "jaccard")
  )
  
  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_cosine_plot",
    clonalOverlap(
      combined, 
      chain = "TRB",
      method = "cosine")
  )
  
  expect_doppelganger( # warning from testthat: Removed 36 rows containing missing values (`geom_text()`).
    "clonalOverlap_coverlap_plot",
    clonalOverlap(
      combined, 
      method = "overlap")
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