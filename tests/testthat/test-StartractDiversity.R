# test script for StartracDiversity.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)
test_obj$Patient <- substr(test_obj$orig.ident,1,3)
test_obj$Type <- substr(test_obj$orig.ident,4,4)


test_that("StartracDiversity works", {
  
  #Ridiculous ggplot warning system can't suppress any internal warnings about NAs
  expect_doppelganger(
    "StartracDiversity_plot",
    suppressWarnings(print(StartracDiversity(test_obj, 
                                             type = "Type", 
                                             group.by = "Patient")))
  )
  
  expect_equal(
    StartracDiversity(test_obj, 
                      type = "Type", 
                      group.by = "Patient",
                      exportTable = TRUE),
    getdata("seuratFunctions", "StartracDiversity_exportTable")
  )
})
