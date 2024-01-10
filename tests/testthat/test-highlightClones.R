# test script for highlightClones.R - testcases are NOT comprehensive!

test_that("highlightClones works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)

  test_obj <- highlightClones(test_obj, 
                              cloneCall = "aa",
                              sequence = "CVVSDNTGGFKTIF_CASSVRRERANTGELFF")
  expect_equal(
    test_obj@meta.data[1:10, 7:14],
    getdata("seuratFunctions", "highlightClones_meta")
  )
})
