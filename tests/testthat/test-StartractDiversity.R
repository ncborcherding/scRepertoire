# test script for StartracDiversity.R - testcases are NOT comprehensive!

test_that("StartracDiversity works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  
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
    getdata("seuratFunctions", "StartracDiversity_exportTable"),
    tolerance = 1e-6
  )
})
