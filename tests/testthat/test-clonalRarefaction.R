# test script for clonalRarefaction.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("clonalRarefaction works", {
  expect_doppelganger(
    "clonalclonalRarefaction_h0_p1_plot",
    clonalRarefaction(combined[1:2], 
                      plot.type = 1,
                      hill.numbers = 0,
                      n.boots = 1)
  )
  
  trial1 <- clonalRarefaction(combined[1:2], 
                              plot.type = 1,
                              hill.numbers = 0,
                              n.boots = 1, 
                              exportTable = TRUE)
  expect_identical(trial1, getdata("visualizations", "clonalRarefaction_h0_p1_exportTable"))
  
  expect_doppelganger(
    "clonalclonalRarefaction_h1_p2_plot",
    clonalRarefaction(combined[3:4], 
                      plot.type = 2,
                      hill.numbers = 1,
                      n.boots = 1)
  )
  
  trial2 <- clonalRarefaction(combined[3:4], 
                              plot.type = 2,
                              hill.numbers = 1,
                              n.boots = 1, 
                              exportTable = TRUE)
  expect_identical(trial2, getdata("visualizations", "clonalRarefaction_h1_p2_exportTable"))
  
  expect_doppelganger(
    "clonalclonalRarefaction_h2_p3_plot",
    clonalRarefaction(combined[5:6], 
                      plot.type = 3,
                      hill.numbers = 2,
                      n.boots = 1)
  )
  
  
  trial3 <- clonalRarefaction(combined[5:6], 
                               plot.type = 3,
                               hill.numbers = 2,
                               n.boots = 1, 
                               exportTable = TRUE)
  expect_identical(trial3, getdata("visualizations", "clonalRarefaction_h2_p3_exportTable"))
  

})