# test script for percentVJ.R - testcases are NOT comprehensive!

test_that("percentVJ works", {

  combined <- getCombined()

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