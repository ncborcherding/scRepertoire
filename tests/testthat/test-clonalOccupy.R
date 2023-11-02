# test script for clonalOccupy.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)


test_that("clonalOccupy works", {

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
