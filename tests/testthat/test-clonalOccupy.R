# test script for clonalOccupy.R - testcases are NOT comprehensive!

test_that("clonalOccupy works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)

  expect_doppelganger(
    "clonalOccupy_default_plot",
    clonalOccupy(test_obj, 
                 x.axis = "seurat_clusters")
  )
  expect_doppelganger(
    "clonalOccupy_proportion_plot",
    clonalOccupy(test_obj, 
                 x.axis = "ident", 
                 proportion = TRUE, 
                 label = FALSE)
  )
  expect_equal(
    clonalOccupy(test_obj, 
                 x.axis = "seurat_clusters", 
                 exportTable = TRUE),
    getdata("seuratFunctions", "clonalOccupy_exportTable")
  )
})
