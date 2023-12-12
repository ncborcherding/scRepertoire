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
  expect_equal(
  clonalCluster(combined[1:2], 
                chain = "TRB", 
                sequence = "aa"),
  getdata("clustering", "clonalCluster_2sample_data")
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
  
  
  BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
  combined.BCR <- combineBCR(BCR, 
                             samples = "S1")
  expect_equal(
      clonalCluster(combined.BCR, 
                    chain = "IGH", 
                    sequence = "aa"),
      getdata("clustering", "clonalCluster_IGHaa_data")
  )
  
})
#TODO Add exportgraph test