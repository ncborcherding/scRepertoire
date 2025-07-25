# test script for StartracDiversity.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Input validation and error handling work correctly", {
  expect_error(
    StartracDiversity(scRep_example, index = c("migr", "tran"), pairwise = "Type"),
    "Pairwise analysis can only be performed for a single index"
  )
  
  expect_error(
    StartracDiversity(scRep_example, index = "expa", pairwise = "Type"),
    "Pairwise analysis is only supported for 'migr' or 'tran' indices."
  )
})


test_that("Output format and structure are correct", {
  # Returns a ggplot object by default
  plot_output <- StartracDiversity(scRep_example, 
                                   type = "Type", 
                                   group.by = "Patient")
  expect_s3_class(plot_output, "ggplot")
  
  # Returns a data.frame when exportTable = TRUE
  table_output <- StartracDiversity(scRep_example, 
                                    type = "Type", 
                                    group.by = "Patient", 
                                    exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Check standard output columns
  expect_true(all(c("group", "cluster", "expa", "migr", "tran") %in% names(table_output)))
})


test_that("Helper functions for entropy calculate correctly", {
  # Test matrix
  mat <- matrix(c(1, 1, 0, 4, 2, 0, 3, 0, 0), nrow = 3, byrow = TRUE)
  
  expected_row_entropy <- c(1, 0.9182958, 0)
  expect_equal(.mrowEntropy(mat), expected_row_entropy, tolerance = 1e-6)

  expected_col_entropy <- c(1.405639, 0.9182958, NA)
  expect_equal(.mcolEntropy(mat), expected_col_entropy, tolerance = 1e-2)
})


test_that("Standard index calculations are mathematically correct", {
  results <- StartracDiversity(scRep_example, 
                               type = "Type", 
                               group.by = "Patient", 
                               exportTable = TRUE)
  

  expect_equal(results$expa[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$expa[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$expa[results$cluster == "3"][1], 0, tolerance = 1e-4)
  
  expect_equal(results$migr[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$migr[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$migr[results$cluster == "3"][1], 0, tolerance = 1e-4)
  
  expect_equal(results$tran[results$cluster == "1"][1], 0, tolerance = 1e-4)
  expect_equal(results$tran[results$cluster == "2"][1], 0, tolerance = 1e-4)
  expect_equal(results$tran[results$cluster == "3"][1], 0.119958, tolerance = 1e-4)
})


