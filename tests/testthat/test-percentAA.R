# test script for percentAA.R - testcases are NOT comprehensive!

test_that("percentAA works", {

  combined <- getCombined()

  set.seed(42)
  expect_doppelganger(
    "percentAA_order_plot",
    percentAA(combined, 
              chain = "TRB", 
              aa.length = 20,
              order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
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