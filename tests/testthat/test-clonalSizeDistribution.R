# test script for clonalSizeDistribution.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                      samples = c("P17B", "P17L", "P18B", "P18L",
                                  "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))



test_that("Output formats (distance matrix vs. plot) are correct", {
  # Test that exportTable = TRUE returns a matrix
  dist_matrix <- clonalSizeDistribution(combined[1:2], exportTable = TRUE)
  expect_true(is.matrix(dist_matrix))
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalSizeDistribution(combined[1:2])
  expect_s3_class(plot_output, "ggplot")
})

test_that("Distance matrix has correct properties", {
  dist_matrix <- clonalSizeDistribution(combined[1:3], exportTable = TRUE)
  num_samples <- length(combined[1:3])
  
  # Matrix should be square and have correct dimensions
  expect_equal(nrow(dist_matrix), num_samples)
  expect_equal(ncol(dist_matrix), num_samples)
  
  # Should have correct row and column names
  expect_equal(rownames(dist_matrix), names(combined)[1:3])
  expect_equal(colnames(dist_matrix), names(combined)[1:3])
  
  # Diagonal should be zero
  expect_true(all(diag(dist_matrix) == 0))
  
  # Matrix should be symmetric
  expect_equal(dist_matrix, t(dist_matrix))
  
  # All values should be non-negative
  expect_true(all(dist_matrix >= 0))
})

test_that("`group.by` parameter correctly aggregates data", {
  dist_matrix_grouped <- clonalSizeDistribution(combined, group.by = "Type", exportTable = TRUE)
  num_groups <- 2
  expect_equal(nrow(dist_matrix_grouped), num_groups)
  expect_equal(ncol(dist_matrix_grouped), num_groups)
  expect_true(all(rownames(dist_matrix_grouped) %in% c("B", "L")))
})

test_that("`method` parameter is passed to hclust correctly", {
  expect_silent(clonalSizeDistribution(combined[1:2], method = "ward.D2"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "complete"))
  expect_silent(clonalSizeDistribution(combined[1:2], method = "single"))
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  expect_silent(
    res_aa <- clonalSizeDistribution(combined[1:2], cloneCall = "aa", exportTable = TRUE)
  )
  expect_true(is.matrix(res_aa))
  
  expect_silent(
    res_tra <- clonalSizeDistribution(combined[1:2], chain = "TRA", exportTable = TRUE)
  )
  expect_true(is.matrix(res_tra))
})

test_that("Plotting output has the correct dendrogram components", {
  plot_output <- clonalSizeDistribution(combined[5:6])
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomSegment")
  expect_s3_class(plot_output$layers[[2]]$geom, "GeomText")
  expect_s3_class(plot_output$layers[[3]]$geom, "GeomPoint")
  expect_s3_class(plot_output$coordinates, "CoordFlip")
  expect_true(plot_output$scales$scales[[1]]$trans$name == "reverse")
})
