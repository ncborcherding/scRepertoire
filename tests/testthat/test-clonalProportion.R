# test script for clonalProportion.R - testcases are NOT comprehensive!

test_that("clonalProportion works", {

  combined <- getCombined()

  expect_doppelganger(
    "clonalProportion_plot",
    clonalProportion(combined, 
                    cloneCall = "gene")
  )
  expect_equal(
    clonalProportion(combined, 
                    cloneCall = "gene", 
                    exportTable = TRUE),
      getdata("visualizations", "clonalProportion_exportTable")
  )
})