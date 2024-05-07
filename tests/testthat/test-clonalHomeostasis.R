# test script for clonalHomeostasis.R - testcases are NOT comprehensive!

test_that("clonalDiversity works", {

  combined <- getCombined()

  expect_doppelganger(
    "clonalHomeostasis_plot",
    clonalHomeostasis(combined, 
                    cloneCall = "gene")
  )
  
  expect_doppelganger(
    "clonalHomeostasis_order_plot",
    clonalHomeostasis(combined, 
                      cloneCall = "gene",
                      order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  expect_equal(
    clonalHomeostasis(combined, 
                      cloneCall = "gene", 
                      exportTable = TRUE),
      getdata("visualizations", "clonalHomeostasis_exportTable")
  )
})