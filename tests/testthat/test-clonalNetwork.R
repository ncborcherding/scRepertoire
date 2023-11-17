# test script for clonalNetwork.R - testcases are NOT comprehensive!

# library(ggraph) # Qile: I believe libraries shouldn't need to be loaded in scripts? Correct me if im wrong.

test_that("clonalNetwork works", {
  library(ggraph)
  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  
  set.seed(42)
  expect_doppelganger( # warning from testthat: Using the `size` aesthetic in this geom was deprecated in ggplot2 3.4.0. Please use `linewidth` in the `default_aes` field and elsewhere instead.
    "clonalNetwork_plot",
    clonalNetwork(test_obj, 
                  reduction = "umap", 
                  group.by = "seurat_clusters", 
                  filter.identity = 3,
                  cloneCall = "aa")
  )
  expect_equal(
    clonalNetwork(test_obj, 
                  reduction = "umap", 
                  group.by = "seurat_clusters", 
                  cloneCall = "aa", 
                  exportTable = TRUE),
    getdata("seuratFunctions", "clonalNetwork_exportTable")
  )
  
  expect_equal(
    clonalNetwork(test_obj, 
                  reduction = "umap", 
                  group.by = "seurat_clusters", 
                  cloneCall = "aa", 
                  exportClones = TRUE),
    getdata("seuratFunctions", "clonalNetwork_exportClones")
  )
  
  
})
