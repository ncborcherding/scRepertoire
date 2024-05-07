# test script for clonalProportion.R - testcases are NOT comprehensive!

test_that("clonalProportion works", {

  combined <- getCombined()

  expect_doppelganger(
    "clonalProportion_plot",
    clonalProportion(combined, 
                    cloneCall = "gene")
  )
  
  
  expect_doppelganger(
    "clonalProportion_order_plot",
    clonalProportion(combined, 
                     cloneCall = "gene",
                     order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  expect_equal(
    clonalProportion(combined, 
                    cloneCall = "gene", 
                    exportTable = TRUE),
      getdata("visualizations", "clonalProportion_exportTable")
  )
})