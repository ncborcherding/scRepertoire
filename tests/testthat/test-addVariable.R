# test script for addVariable.R - testcases are NOT comprehensive!

test_that("addVariable works", {
  expect_equal(addVariable(getCombined(), 
                           variable.name = "Type", 
                           variables = rep(c("B", "L"), 4)),
               getdata("processing", "addVariable_data"))
})