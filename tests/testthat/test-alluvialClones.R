# test script for alluvialClones.R - testcases are NOT comprehensive!

test_that("alluvialClones works", {

  data("scRep_example")
  test_obj <- combineExpression(getdata("combineContigs", "combined"), scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  
  expect_doppelganger(
    "alluvialClones_plot",
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident"), 
                   color = "Type")
  )
  
  expect_doppelganger(
    "alluvialClones_nocolor_plot",
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident"), 
                   color = NULL)
  )
  
  expect_doppelganger(
    "alluvialClones_facet_plot",
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident"), 
                   color = "ident", 
                   facet = "Type")
  )
  
  
  expect_equal(
    alluvialClones(test_obj, 
                   cloneCall = "aa", 
                   y.axes = c("Patient", "ident", "Type"), 
                   exportTable = TRUE),
    getdata("seuratFunctions", "alluvialClones_exportTable")
  )
})
