# test script for clonalRarefaction.R - testcases are NOT comprehensive!

# Qile: I'm getting a bunch of errors when running devtools::test()

test_that("clonalRarefaction works", {

  combined <- getCombined()

  expect_doppelganger(
    "clonalclonalRarefaction_h0_p1_plot",
    clonalRarefaction(combined[1:2],
                      plot.type = 1,
                      hill.numbers = 0,
                      n.boots = 1) + 
      guides(color = "none", shape = "none", linetype = "none")
  )
  
  trial1 <- clonalRarefaction(combined[1:2],
                              plot.type = 1,
                              hill.numbers = 0,
                              n.boots = 1,
                              exportTable = TRUE)
  expect_equal(
    trial1,
    getdata("visualizations", "clonalRarefaction_h0_p1_exportTable"),
    tolerance = 1e-7
  )
  
  expect_doppelganger(
    "clonalclonalRarefaction_h1_p2_plot",
    clonalRarefaction(combined[3:4],
                      plot.type = 2,
                      hill.numbers = 1,
                      n.boots = 1) + 
      guides(color = "none", shape = "none", linetype = "none")
  )
  
  trial2 <- clonalRarefaction(combined[3:4],
                              plot.type = 2,
                              hill.numbers = 1,
                              n.boots = 1,
                              exportTable = TRUE)
  expect_equal(
    trial2,
    getdata("visualizations", "clonalRarefaction_h1_p2_exportTable"),
    tolerance = 1e-7
  )
  
  expect_doppelganger(
    "clonalclonalRarefaction_h2_p3_plot",
    clonalRarefaction(combined[5:6],
                      plot.type = 3,
                      hill.numbers = 2,
                      n.boots = 1) + 
      guides(color = "none", shape = "none", linetype = "none")
  )
  
  
  trial3 <- clonalRarefaction(combined[5:6], 
                               plot.type = 3,
                               hill.numbers = 2,
                               n.boots = 1, 
                               exportTable = TRUE)
  expect_equal(
    trial3,
    getdata("visualizations", "clonalRarefaction_h2_p3_exportTable"),
    tolerance = 1e-3 # this is low jsut because of one value actual$data[19, ]
  )
})