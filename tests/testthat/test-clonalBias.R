# test script for clonalBias.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)
test_obj$Patient <- substr(test_obj$orig.ident,1,3)
test_obj$Type <- substr(test_obj$orig.ident,4,4)


test_that("clonalBias works", {
  
  set.seed(42)
  expect_doppelganger(
    "clonalBias_plot",
    clonalBias(test_obj, 
               cloneCall = "aa", 
               split.by = "Patient", 
               group.by = "seurat_clusters",
               n.boots = 1, 
               min.expand = 2)
  )
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
