# test script for StartracDiversity.R - testcases are NOT comprehensive!

test_that("StartracDiversity works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  
  expect_doppelganger(
    "StartracDiversity_plot",
    StartracDiversity(test_obj, 
                      type = "Type", 
                       group.by = "Patient")
  )
})
