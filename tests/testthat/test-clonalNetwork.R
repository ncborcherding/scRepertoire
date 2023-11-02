# test script for clonalNetwork.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
library(ggraph)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)



test_that("clonalNetwork works", {
  
  set.seed(42)
  expect_doppelganger(
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
