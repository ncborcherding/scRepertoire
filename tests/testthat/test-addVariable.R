# test script for addVariable.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))

test_that("Function correctly adds a new variable column", {
  variable_values <- rep(c("B_cell", "L_cell"), 4)
  
  combined_new <- addVariable(combined, 
                              variable.name = "Type", 
                              variables = variable_values)
  
  expect_type(combined_new, "list")
  expect_length(combined_new, length(combined))
  
  for (i in seq_along(combined_new)) {
    df <- combined_new[[i]]
    expect_true("Type" %in% names(df))
    expect_equal(unique(df$Type), variable_values[i])
  }
})

test_that("Function throws an error with mismatched lengths", {
  mismatched_variables <- c("A", "B", "C")
  error_message <- "Make sure the variables match the length of the contig list"
  expect_error(addVariable(combined, 
                           variable.name = "Mismatch", 
                           variables = mismatched_variables),
               error_message, fixed = TRUE)
})

test_that("Function can handle different data types", {
  numeric_vars <- 1:length(combined)
  combined_numeric <- addVariable(combined, 
                                  variable.name = "Batch", 
                                  variables = numeric_vars)
  expect_true("Batch" %in% names(combined_numeric[[1]]))
  expect_type(combined_numeric[[1]]$Batch, "integer")
  expect_equal(unique(combined_numeric[[3]]$Batch), 3)
  factor_vars <- factor(rep(c("Group1", "Group2"), 4))
  combined_factor <- addVariable(combined, 
                                 variable.name = "Group", 
                                 variables = factor_vars)
  
  expect_true("Group" %in% names(combined_factor[[1]]))
  expect_s3_class(combined_factor[[1]]$Group, "factor")
})

test_that("Function can overwrite an existing column", {
  combined_initial <- addVariable(combined, 
                                  variable.name = "Status", 
                                  variables = rep("Original", length(combined)))
  expect_equal(unique(combined_initial[[1]]$Status), "Original")
  combined_overwritten <- addVariable(combined_initial, 
                                      variable.name = "Status", 
                                      variables = rep("Updated", length(combined)))
  for (df in combined_overwritten) {
    expect_equal(unique(df$Status), "Updated")
  }
})

test_that("Function works with a single-element list", {
  single_element_list <- list(P17B = combined$P17B)
  single_variable <- "SingleSample"
  result <- addVariable(single_element_list, 
                        variable.name = "Experiment", 
                        variables = single_variable)
  
  expect_length(result, 1)
  expect_true("Experiment" %in% names(result$P17B))
  expect_equal(unique(result$P17B$Experiment), "SingleSample")
})
