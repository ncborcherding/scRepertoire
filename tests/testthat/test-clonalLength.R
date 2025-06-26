# test script for clonalLength.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Input validation for `cloneCall` works", {
  error_message <- "Please make a selection of the type of CDR3 sequence to analyze by using `cloneCall`"
  expect_error(clonalLength(combined, cloneCall = "gene"), error_message, fixed = TRUE)
  expect_error(clonalLength(combined, cloneCall = "strict"), error_message, fixed = TRUE)
})

test_that("Output formats (table vs. plot) are correct", {
  table_output <- clonalLength(combined, exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  expect_true(all(c("length", "values") %in% names(table_output)))
  plot_output <- clonalLength(combined)
  expect_s3_class(plot_output, "ggplot")
})

test_that("`scale` parameter correctly switches plot type", {
  # Test with scale = FALSE (default) for a bar plot
  plot_bar <- clonalLength(combined, scale = FALSE)
  expect_s3_class(plot_bar$layers[[1]]$geom, "GeomBar")
  expect_equal(plot_bar$labels$y, "Number of CDR3 (AA)")
  
  # Test with scale = TRUE for a density plot
  plot_density <- clonalLength(combined, scale = TRUE)
  expect_s3_class(plot_density$layers[[1]]$geom, "GeomDensity")
  expect_equal(plot_density$labels$y, "Density of CDR3 (AA)")
})

test_that("`group.by` parameter functions correctly", {
  # Test exporting a table grouped by "Type"
  table_grouped <- clonalLength(combined, group.by = "Type", exportTable = TRUE)
  expect_true("Type" %in% names(table_grouped))
  expect_true("values" %in% names(table_grouped))
  expect_true(all(unique(table_grouped$Type) %in% c("B", "L")))
  
  # Test plotting grouped by "Type"
  plot_grouped <- clonalLength(combined, group.by = "Type")
  expect_s3_class(plot_grouped, "ggplot")
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  # Test with cloneCall = "nt"
  plot_nt <- clonalLength(combined, cloneCall = "nt")
  expect_s3_class(plot_nt, "ggplot")
  expect_equal(plot_nt$labels$y, "Number of CDR3 (NT)")
  table_nt <- clonalLength(combined, cloneCall = "nt", exportTable = TRUE)
  table_aa <- clonalLength(combined, cloneCall = "aa", exportTable = TRUE)
  expect_gt(mean(table_nt$length), mean(table_aa$length))
  
  # Test single chain selection
  expect_silent(plot_tra <- clonalLength(combined, chain = "TRA"))
  expect_s3_class(plot_tra, "ggplot")
  table_tra <- clonalLength(combined, chain = "TRA", exportTable = TRUE)
  expect_s3_class(table_tra, "data.frame")
  expect_lt(nrow(table_tra), nrow(table_aa))
})
