# test script for clonalLength.R - testcases are NOT comprehensive!

test_that("clonalLength works", {
  expect_doppelganger(
    "clonalLength_both_chain_plot", 
    clonalLength(getCombined(), chain = "both") 
  )
  
  expect_doppelganger(
    "clonalLength_both_chain_order_plot", 
    clonalLength(getCombined(), 
                 chain = "both",
                 order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L")) 
  )
  
  expect_doppelganger(
    "clonalLength_scaled_plot", 
    clonalLength(getCombined(), 
                 scale = TRUE) 
  )
  
  expect_doppelganger(
    "clonalLength_TRB_plot", 
    clonalLength(getCombined(), 
                 chain = "TRB") 
  )
  
  expect_doppelganger(
    "clonalLength_TRA_plot", 
    clonalLength(getCombined(), 
                 chain = "TRA") 
  )
  
  combined <- addVariable(getCombined(), 
                          variable.name = "Type", 
                          variables = rep(c("B", "L"), 4))
  expect_doppelganger(
    "clonalLength_groupby_plot",
    clonalLength(combined, 
                    cloneCall = "nt",
                    group.by = "Type")
  )
})

