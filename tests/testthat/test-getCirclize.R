# test script for getCirclize.R - testcases are NOT comprehensive!

test_that("getCirclize works", {
    data("scRep_example")
    test_obj <- combineExpression(getdata("combineContigs", "combined"), scRep_example)
    test_obj$Patient <- substr(test_obj$orig.ident,1,3)

    expect_equal(
      getCirclize(test_obj, 
                  cloneCall = "aa"),
      getdata("seuratFunctions", "getCirclize_default_output")
    )    
    
    expect_equal(
      getCirclize(test_obj, 
                  cloneCall = "nt",
                  proportion = TRUE,
                  group.by = "Patient"),
      getdata("seuratFunctions", "getCirclize_proportion_output")
    )    
})