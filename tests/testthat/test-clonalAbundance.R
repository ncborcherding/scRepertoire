# test script for clonalAbundance.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))

test_that("exportTable = TRUE returns a data.frame", {
  table_output <- clonalAbundance(combined, exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  expect_true(all(c("Abundance", "values") %in% names(table_output)))
})

test_that("Default behavior returns a ggplot object", {
  plot_output <- clonalAbundance(combined)
  expect_s3_class(plot_output, "ggplot")
})

test_that("scale = TRUE produces a density plot", {
  plot_output <- clonalAbundance(combined, scale = TRUE)
  expect_s3_class(plot_output, "ggplot")
  expect_equal(plot_output$labels$y, "Density of Clones")
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomDensity")
})

test_that("scale = FALSE produces a line plot", {
  plot_output <- clonalAbundance(combined, scale = FALSE)
  expect_s3_class(plot_output, "ggplot")
  expect_equal(plot_output$labels$y, "Number of Clones")
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomLine")
})

test_that("group.by correctly groups the data", {
  table_output_grouped <- clonalAbundance(combined, group.by = "Type", exportTable = TRUE)
  expect_s3_class(table_output_grouped, "data.frame")
  expect_true("Type" %in% names(table_output_grouped))
  expect_lt(length(unique(table_output_grouped$patient)), length(combined))
  
  plot_output_grouped <- clonalAbundance(combined, group.by = "Type")
  expect_s3_class(plot_output_grouped, "ggplot")
  expect_equal(plot_output_grouped$labels$colour, "Type")
})

test_that("cloneCall handles different arguments", {
  # Test with cloneCall = "gene"
  expect_s3_class(clonalAbundance(combined, cloneCall = "gene", exportTable = TRUE), "data.frame")
  expect_s3_class(clonalAbundance(combined, cloneCall = "gene"), "ggplot")
  
  # Test with cloneCall = "nt"
  expect_s3_class(clonalAbundance(combined, cloneCall = "nt", exportTable = TRUE), "data.frame")
  expect_s3_class(clonalAbundance(combined, cloneCall = "nt"), "ggplot")
  
  # Test with cloneCall = "aa"
  expect_s3_class(clonalAbundance(combined, cloneCall = "aa", exportTable = TRUE), "data.frame")
  expect_s3_class(clonalAbundance(combined, cloneCall = "aa"), "ggplot")
})

test_that("chain selection works correctly", {
  table_TRA <- clonalAbundance(combined, chain = "TRA", exportTable = TRUE)
  expect_s3_class(table_TRA, "data.frame")
  
  # Test with chain = "TRB"
  table_TRB <- clonalAbundance(combined, chain = "TRB", exportTable = TRUE)
  expect_s3_class(table_TRB, "data.frame")
  
  # Test with chain = "both" (default)
  table_both <- clonalAbundance(combined, chain = "both", exportTable = TRUE)
  expect_s3_class(table_both, "data.frame")
  
  expect_lte(nrow(table_TRA), nrow(table_both))
  expect_lte(nrow(table_TRB), nrow(table_both))
})
