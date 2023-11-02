# test script for alluvialClones.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)
test_obj$Patient <- substr(test_obj$orig.ident,1,3)
test_obj$Type <- substr(test_obj$orig.ident,4,4)


test_that("alluvialClones works", {
  
  expect_doppelganger(
    "alluvialClones_plot",
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident"), 
                   color = "Type")
  )
  
  expect_equal(
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident", "Type"), 
                   exportTable = TRUE),
    getdata("seuratFunctions", "alluvialClones_exportTable")
  )
})
