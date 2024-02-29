# test script for clonalOverlay.R - testcases are NOT comprehensive!

test_that("clonalOverlay works", {

  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  test_obj$Patient <- substr(test_obj$orig.ident, 1,3)
  
  set.seed(42)
  expect_doppelganger(
    "clonalOverlay_plot",
    clonalOverlay(test_obj, 
                 reduction = "umap",
                 cutpoint = 1, 
                 bins = 10, 
                 facet.by = "Patient") + 
      guides(color = "none")
  )
  
  set.seed(42)
  expect_doppelganger(
    "clonalOverlay_clonalProportion_plot",
    clonalOverlay(test_obj, 
                  reduction = "umap",
                  cut.category = "clonalProportion",
                  cutpoint = 0.2, 
                  bins = 10, 
                  facet.by = "Patient") + 
      guides(color = "none")
  )
  
})
