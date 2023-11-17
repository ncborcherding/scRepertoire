# test script for clonalCompare.R - testcases are NOT comprehensive!

test_that("clonalCompare works", {
  
  combined <- getCombined()

  expect_doppelganger(
    "clonalCompare_alluvial_plot",
    clonalCompare(
      combined, 
      top.clones = 10, 
      samples = c("P17B", "P17L"), 
      cloneCall="aa", 
      graph = "alluvial"
    )
  )
  
  expect_doppelganger(
    "clonalCompare_area_plot",
    clonalCompare(
      combined, 
      top.clones  = 10, 
      samples = c("P17B", "P17L"), 
      cloneCall="aa", 
      graph = "area")
  )
})
