# test script for alluvialClones.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Input validation works correctly", {
  expect_error(alluvialClones(scRep_example, y.axes = NULL),
               "Make sure you have selected the variable\\(s\\) to visualize")
})

test_that("Output formats (table vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame in lodes format
  table_output <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"), exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  expect_true(all(c("x", "stratum", "alluvium") %in% names(table_output)))
  
  # Test that the default behavior returns a ggplot object
  plot_output <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"))
  expect_s3_class(plot_output, "ggplot")
})

test_that("Aesthetic parameters (`color`, `alpha`, `facet`) are handled correctly", {
  # Test with color as a column
  plot_color <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"), color = "Type")
  expect_equal(plot_color$labels$fill, "Type")
  
  # Test with color as a specific clone
  plot_color_clone <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"), color = "cloneA")
  expect_equal(plot_color_clone$labels$fill, "clone(s)")
  
  # Test with facet
  plot_facet <- alluvialClones(scRep_example, y.axes = c("ident", "Type"), facet = "Patient")
  expect_s3_class(plot_facet$facet, "FacetWrap")
})

test_that("Lodes data frame is generated correctly", {
  lodes_df <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"), exportTable = TRUE)
  axis_labels <- unique(as.character(lodes_df$x))
  expect_equal(axis_labels, c("Patient", "Type"))
})

test_that("Plot structure has correct geoms", {
  # No color or alpha
  plot_plain <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"))
  expect_true("GeomStratum" %in% sapply(plot_plain$layers, function(l) class(l$geom)[1]))
  expect_true("GeomAlluvium" %in% sapply(plot_plain$layers, function(l) class(l$geom)[1]))
  
  # With color
  plot_color <- alluvialClones(scRep_example, y.axes = c("Patient", "Type"), color = "Type")
  expect_true("GeomFlow" %in% sapply(plot_color$layers, function(l) class(l$geom)[1]))
  
  # All plots should have geom_text
  expect_true("GeomText" %in% sapply(plot_plain$layers, function(l) class(l$geom)[1]))
})
