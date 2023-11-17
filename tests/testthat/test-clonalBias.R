# test script for clonalBias.R - testcases are NOT comprehensive!

test_that("clonalBias works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  
  #TODO clonalBias plot check
  set.seed(42)
  expect_equal(
    clonalBias(test_obj, 
               cloneCall = "aa", 
               split.by = "Patient", 
               group.by = "seurat_clusters",
               n.boots = 1, 
               min.expand = 2, 
              exportTable = TRUE),
    getdata("seuratFunctions", "clonalBias_exportTable")
  )
})
