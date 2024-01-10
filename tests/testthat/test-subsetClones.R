# test script for subsetClones.R - testcases are NOT comprehensive!

test_that("subsetClones works", {
  expect_equal(subsetClones(getCombined(), name = "sample", variables = c("P17B")),
               getdata("processing", "subsetClones_data"))
})