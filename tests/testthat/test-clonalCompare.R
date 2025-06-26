# test script for clonalCompare.R - testcases are NOT comprehensive!


# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))


test_that("exportTable = TRUE returns a data.frame", {
  table_output <- clonalCompare(combined, top.clones = 5, exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  expect_true(all(c("clones", "Proportion", "Sample") %in% names(table_output)))
})

test_that("Default behavior returns a ggplot object (alluvial)", {
  plot_output <- clonalCompare(combined, top.clones = 5)
  expect_s3_class(plot_output, "ggplot")
  geoms <- sapply(plot_output$layers, function(l) class(l$geom)[1])
  expect_true("GeomStratum" %in% geoms)
  expect_true("GeomFlow" %in% geoms)
})

test_that("graph = 'area' produces an area plot", {
  plot_output <- clonalCompare(combined, top.clones = 5, graph = "area")
  expect_s3_class(plot_output, "ggplot")
  geoms <- sapply(plot_output$layers, function(l) class(l$geom)[1])
  expect_true("GeomArea" %in% geoms)
})

test_that("proportion = FALSE uses 'Count' and proportion = TRUE uses 'Proportion'", {
  # Test with proportion = FALSE (raw counts)
  table_count <- clonalCompare(combined, top.clones = 5, proportion = FALSE, exportTable = TRUE)
  expect_true("Count" %in% names(table_count))
  expect_false("Proportion" %in% names(table_count))
  expect_true(all(table_count$Count >= 1))
  expect_type(table_count$Count, "integer")
  
  # Test with proportion = TRUE (relative proportions)
  table_prop <- clonalCompare(combined, top.clones = 5, proportion = TRUE, exportTable = TRUE)
  expect_true("Proportion" %in% names(table_prop))
  expect_false("Count" %in% names(table_prop))
  expect_true(all(table_prop$Proportion > 0 & table_prop$Proportion <= 1))
})

test_that("Filtering by `samples`, `clones`, and `top.clones` works", {
  # Get a list of top clones to use for the `clones` parameter test
  all_clones <- clonalCompare(combined, exportTable = TRUE, top.clones = 2)
  clones_to_filter <- unique(all_clones$clones)
  
  # Test filtering by specific clones
  table_clones <- clonalCompare(combined, clones = as.character(clones_to_filter), exportTable = TRUE)
  expect_true(all(unique(table_clones$clones) %in% clones_to_filter))
  
  # Test filtering by specific samples
  samples_to_filter <- c("P17B", "P18L")
  table_samples <- clonalCompare(combined, top.clones = 5, samples = samples_to_filter, exportTable = TRUE)
  expect_true(all(unique(table_samples$Sample) %in% samples_to_filter))
  
  # Test that `clones` takes precedence over `top.clones`
  table_precedence <- clonalCompare(combined, clones = as.character(clones_to_filter), top.clones = 10, exportTable = TRUE)
  expect_true(all(unique(table_precedence$clones) %in% clones_to_filter))
  expect_length(unique(table_precedence$clones), length(clones_to_filter))
  
  # Test error when filtering results in no clones
  expect_error(clonalCompare(combined, clones = c("nonexistent_clone_1", "nonexistent_clone_2")),
               regexp = "Please reasses the filtering strategies here")
  
})

test_that("group.by correctly groups the data", {
  table_grouped <- clonalCompare(combined, top.clones = 3, group.by = "Type", exportTable = TRUE)
  expect_s3_class(table_grouped, "data.frame")
  expect_true(all(unique(table_grouped$Sample) %in% c("B", "L")))
})

test_that("relabel.clones = TRUE renames clones numerically", {
  table_relabel <- clonalCompare(combined, top.clones = 4, relabel.clones = TRUE, exportTable = TRUE)
  expect_true(all(str_detect(unique(table_relabel$clones), "^Clone: \\d+$")))
  expect_true("original.clones" %in% names(table_relabel))
})

test_that("highlight.clones modifies colors correctly", {
  # Get top clones to use for highlighting
  clones_to_highlight <- clonalCompare(combined, top.clones = 2, exportTable = TRUE)$clones
  plot_highlight <- clonalCompare(combined, top.clones = 10, highlight.clones = clones_to_highlight)
  
  # Build the plot to inspect its components
  plot_build <- ggplot_build(plot_highlight)
  plot_data <- plot_build$data[[1]] # First layer data
  
  # Check the number of unique fill colors
  unique_colors <- unique(plot_data$fill)
  expect_length(unique_colors, length(unique(clones_to_highlight)) + 1)
  expect_true("grey" %in% unique_colors)
})
