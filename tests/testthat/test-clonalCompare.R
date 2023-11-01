# test script for clonalCompare.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))


test_that("clonalCompare works", {
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
