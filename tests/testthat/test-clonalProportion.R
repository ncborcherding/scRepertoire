# test script for clonalProportion.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Output formats (matrix vs. plot) are correct", {
  mat_output <- clonalProportion(combined, exportTable = TRUE)
  expect_true(is.matrix(mat_output))
  plot_output <- clonalProportion(combined)
  expect_s3_class(plot_output, "ggplot")
})

test_that("The returned matrix has correct dimensions and properties", {
  mat_output <- clonalProportion(combined, exportTable = TRUE)
  expect_equal(nrow(mat_output), length(combined)) # Rows should equal number of samples
  expect_equal(ncol(mat_output), 6)               # Default has 6 bins
  total_clones_per_sample <- sapply(combined, function(x) nrow(x[!is.na(x$CTstrict),]))
  expect_equal(rowSums(mat_output), total_clones_per_sample, tolerance = 1e-9)
  expected_colnames <- c("[1:10]", "[11:100]", "[101:1000]", "[1001:10000]", "[10001:30000]", "[30001:1e+05]")
  expect_equal(colnames(mat_output), expected_colnames)
})

test_that("Custom `clonalSplit` parameter works as expected", {
  # Define custom bins
  custom_bins <- c(5, 50, 500)
  mat_custom <- clonalProportion(combined, clonalSplit = custom_bins, exportTable = TRUE)
  
  # Check dimensions with custom bins
  expect_equal(ncol(mat_custom), 3)
  
  # Check custom column names
  expected_custom_colnames <- c("[1:5]", "[6:50]", "[51:500]")
  expect_equal(colnames(mat_custom), expected_custom_colnames)
})

test_that("`group.by` parameter correctly aggregates the data", {
  mat_grouped <- clonalProportion(combined, group.by = "Type", exportTable = TRUE)
  expect_equal(nrow(mat_grouped), 2)
  expect_true(all(rownames(mat_grouped) %in% c("B", "L")))
  grouped_list <- scRepertoire:::.groupList(combined, "Type")
  total_clones_per_group <- sapply(grouped_list, function(x) nrow(x[!is.na(x$CTstrict),]))
  expect_equal(rowSums(mat_grouped), total_clones_per_group, tolerance = 1e-9)
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  mat_gene <- clonalProportion(combined, cloneCall = "gene", exportTable = TRUE)
  expect_true(is.matrix(mat_gene))
  expect_equal(nrow(mat_gene), length(combined))
  mat_tra <- clonalProportion(combined, chain = "TRA", exportTable = TRUE)
  expect_true(is.matrix(mat_tra))
  expect_equal(nrow(mat_tra), length(combined))
})

test_that("Plotting output has the correct components", {
  # Generate the plot
  plot_output <- clonalProportion(combined)
  
  # Check for the correct geometry and position
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomBar")
  expect_s3_class(plot_output$layers[[1]]$position, "PositionFill")
  
  # Check the labels
  expect_equal(plot_output$labels$y, "Occupied Repertoire Space")
  expect_equal(plot_output$labels$x, "Samples")
  expect_equal(plot_output$labels$fill, "Clonal Indices")
})
