# test script for subsetClones.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list,
                      samples = c("P17B", "P17L", "P18B", "P18L",
                                  "P19B", "P19L", "P20B", "P20L"))


test_that("Subsetting to a single variable works correctly", {
  subset_single <- subsetClones(combined, name = "sample", variables = "P17B")
  expect_type(subset_single, "list")
  expect_length(subset_single, 1)
  expect_equal(names(subset_single), "P17B")
  expect_identical(subset_single[[1]], combined[["P17B"]])
})

test_that("Subsetting to multiple variables works correctly", {
  vars_to_subset <- c("P17L", "P19B", "P20L")
  subset_multi <- subsetClones(combined, name = "sample", variables = vars_to_subset)
  expect_type(subset_multi, "list")
  expect_length(subset_multi, 3)
  expect_equal(names(subset_multi), vars_to_subset)
  expect_identical(subset_multi, combined[vars_to_subset])
})

test_that("Returns an empty list when no variables match", {
  subset_empty <- subsetClones(combined, name = "sample", variables = "nonexistent_sample")
  expect_type(subset_empty, "list")
  expect_length(subset_empty, 0)
  expect_equal(subset_empty, list())
})

test_that("Subsetting works with other added metadata columns", {
  # Add a new variable to the combined object for testing
  combined_with_type <- addVariable(combined, 
                                    variable.name = "Type", 
                                    variables = rep(c("B_cell", "L_cell"), 4))
  subset_type <- subsetClones(combined_with_type, name = "Type", variables = "B_cell")
  expect_length(subset_type, 4)
  expected_names <- c("P17B", "P18B", "P19B", "P20B")
  expect_equal(names(subset_type), expected_names)
})

test_that("Function fails gracefully with an invalid column name", {
  expect_error(subsetClones(combined, name = "invalid_column_name", variables = "any_value"))
})
