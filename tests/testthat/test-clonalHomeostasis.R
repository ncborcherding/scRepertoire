# test script for clonalHomeostasis.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Output formats (matrix vs. plot) are correct", {
  mat_output <- clonalHomeostasis(combined, exportTable = TRUE)
  expect_true(is.matrix(mat_output))
  plot_output <- clonalHomeostasis(combined)
  expect_s3_class(plot_output, "ggplot")
})

test_that("The returned matrix has correct dimensions and properties", {
  mat_output <- clonalHomeostasis(combined, exportTable = TRUE)
  expect_equal(nrow(mat_output), length(combined)) # Rows should equal number of samples
  expect_equal(ncol(mat_output), 5)               # Default has 5 bins
  expect_true(all(abs(rowSums(mat_output) - 1) < 1e-9))
  expected_colnames <- c("Rare (0 < X <= 1e-04)",
                         "Small (1e-04 < X <= 0.001)",
                         "Medium (0.001 < X <= 0.01)",
                         "Large (0.01 < X <= 0.1)",
                         "Hyperexpanded (0.1 < X <= 1)")
  expect_equal(colnames(mat_output), expected_colnames)
})

test_that("Custom `cloneSize` parameter works as expected", {
  custom_bins <- c(Small = 0.01, Medium = 0.2, Large = 1)
  mat_custom <- clonalHomeostasis(combined, cloneSize = custom_bins, exportTable = TRUE)
  expect_equal(ncol(mat_custom), 3)
  expected_custom_colnames <- c("Small (0 < X <= 0.01)",
                                "Medium (0.01 < X <= 0.2)",
                                "Large (0.2 < X <= 1)")
  expect_equal(colnames(mat_custom), expected_custom_colnames)
  expect_true(all(abs(rowSums(mat_custom) - 1) < 1e-9))
})

test_that("`group.by` parameter correctly aggregates the data", {
  mat_grouped <- clonalHomeostasis(combined, group.by = "Type", exportTable = TRUE)
  expect_equal(nrow(mat_grouped), 2)
  expect_true(all(rownames(mat_grouped) %in% c("B", "L")))
  expect_true(all(abs(rowSums(mat_grouped) - 1) < 1e-9))
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  mat_aa <- clonalHomeostasis(combined, cloneCall = "aa", exportTable = TRUE)
  expect_true(is.matrix(mat_aa))
  expect_equal(nrow(mat_aa), length(combined))
  mat_trb <- clonalHomeostasis(combined, chain = "TRB", exportTable = TRUE)
  expect_true(is.matrix(mat_trb))
  expect_equal(nrow(mat_trb), length(combined))
})

test_that("Plotting output has the correct components", {
  plot_output <- clonalHomeostasis(combined)
  
  # Check for the correct geometry and position
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomBar")
  expect_s3_class(plot_output$layers[[1]]$position, "PositionFill")
  
  # Check the labels
  expect_equal(plot_output$labels$y, "Relative Abundance")
  expect_equal(plot_output$labels$x, "Samples")
  expect_equal(plot_output$labels$fill, "Clonal Group")
})
