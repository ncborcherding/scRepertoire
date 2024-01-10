# test script for clonalHomeostasis.R - testcases are NOT comprehensive!

test_that("clonalDiversity works", {

  combined <- getCombined()

  expect_doppelganger(
    "clonalHomeostasis_plot",
    clonalHomeostasis(combined, 
                    cloneCall = "gene")
  )
  expect_equal(
    clonalHomeostasis(combined, 
                      cloneCall = "gene", 
                      exportTable = TRUE),
      getdata("visualizations", "clonalHomeostasis_exportTable")
  )
})