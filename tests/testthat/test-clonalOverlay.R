# test script for clonalOverlay.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)
test_obj$Patient <- substr(test_obj$orig.ident, 1,3)


test_that("clonalOverlay works", {
  
  set.seed(42)
  expect_doppelganger(
    "clonalOverlay_plot",
    clonalOverlay(test_obj, 
                 reduction = "umap",
                 freq.cutpoint = 1, 
                 bins = 10, 
                 facet.by = "Patient")
  )
  
})
