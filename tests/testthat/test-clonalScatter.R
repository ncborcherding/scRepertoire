# test script for clonalScatter.R - testcases are NOT comprehensive!

test_that("clonalScatter works", {
  set.seed(42)
  expect_doppelganger(
    "clonalScatter_proportion_plot",
    clonalScatter(getCombined(), 
                  cloneCall = "gene", 
                  y.axis = "P17B", 
                  x.axis = "P17L")
  )
  
  expect_doppelganger(
    "clonalScatter_raw_plot",
    clonalScatter(getCombined(), 
                  cloneCall = "aa", 
                  y.axis = "P17B", 
                  x.axis = "P17L", 
                  graph = "count")
  )
  
  #exported.graph <- clonalScatter(combined, 
  #                                cloneCall = "aa", 
  #                                chain = "TRB",
  #                                y.axis = "P17B", 
  #                                x.axis = "P17L", 
  #                                exportTable = TRUE)
  #exported.graph <- exported.graph[!grepl("\\;", exported.graph$Var1),]
 # expect_equal(
 #   exported.graph,
 #   getdata("visualizations", "clonalScatter_exportTable"),
  #  tolerance = 1e-4
  #)
  
})
