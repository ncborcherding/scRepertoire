# test script for clonalCluster.R - testcases are NOT comprehensive!

test_that("clonalCluster works", {
 
  data("scRep_example")
  combined <- getCombined()
  test_obj <- combineExpression(combined, scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  
  set.seed(42)
  expect_equal(
    clonalCluster(combined[[1]], 
                  chain = "TRB", 
                  sequence = "aa"),
    getdata("clustering", "clonalCluster_TRBaa_data")
  )
  
  
  set.seed(42)
  test_obj <- clonalCluster(test_obj, 
                      chain = "TRB", 
                      sequence = "aa", 
                      group.by = "Patient")
  expect_equal(
    test_obj@meta.data[, 7:16],
    getdata("clustering", "clonalCluster_TRBaa_metadata")
  )
})
#TODO Add exportgraph test