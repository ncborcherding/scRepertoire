# test script for percentVJ.R - testcases are NOT comprehensive!

test_that("percentVJ works", {

  combined <- getCombined()

  set.seed(42)
  expect_doppelganger(
    "percentVJ_plot",
    percentVJ(combined, 
              chain = "TRB")
  )
  
  expect_doppelganger(
    "percentVJ_order_plot",
    percentVJ(combined, 
              chain = "TRB",
              order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  set.seed(42)
  expect_equal(
    percentVJ(combined, 
              chain = "TRB", 
              exportTable = TRUE),
      getdata("visualizations", "percentVJ_exportTable")
  )
})