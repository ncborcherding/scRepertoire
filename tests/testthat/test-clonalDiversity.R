# test script for clonalDiversity.R - testcases are NOT comprehensive!

test_that("clonalDiversity works", {

  combined <- getCombined()

  set.seed(42)
  expect_doppelganger(
    "clonalDiversity_plot",
    clonalDiversity(combined, 
                    cloneCall = "gene", 
                    n.boots = 1)
  )
  set.seed(42)
  expect_equal(
      clonalDiversity(combined, 
                      cloneCall = "gene", 
                      skip.boots = TRUE, 
                      metrics = c("norm.entropy", "shannon"),
                      exportTable = TRUE),
      getdata("visualizations", "clonalDiversity_exportTable")
  )
})