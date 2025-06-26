# test script for clonalOccupy.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Output formats (table vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame
  table_output <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalOccupy(scRep_example, x.axis = "seurat_clusters")
  expect_s3_class(plot_output, "ggplot")
})

test_that("Data aggregation and counting are correct", {
  table_output <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", exportTable = TRUE)
  c1_single_count <- table_output[table_output$seurat_clusters == "1" & table_output$cloneSize == "Large (0.01 < X <= 0.1)", "n"]
  c1_small_count <- table_output[table_output$seurat_clusters == "1" & table_output$cloneSize == "Medium (0.001 < X <= 0.01)", "n"]
  expect_equal(c1_single_count, 34)
  expect_equal(c1_small_count, 41)
  
})

test_that("`proportion` parameter works correctly", {
  table_prop <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", proportion = TRUE, exportTable = TRUE)
  expect_true(all(c("total", "prop") %in% names(table_prop)))
  total_props <- table_prop %>%
    group_by(seurat_clusters) %>%
    summarise(sum_prop = sum(prop))
  expect_true(all(abs(total_props$sum_prop - 1) < 1e-9))
  plot_prop <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", proportion = TRUE)
  expect_equal(plot_prop$labels$y, "Proportion of Cells")
})

test_that("`facet.by` parameter works correctly", {
  # Test table output with facet
  table_facet <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", facet.by = "Type", exportTable = TRUE)
  expect_true("Type" %in% names(table_facet))
  
  # Test plot output with facet
  plot_facet <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", facet.by = "Type")
  expect_s3_class(plot_facet$facet, "FacetGrid")
})

test_that("`label` parameter correctly adds geom_text", {
  # Test with label = TRUE (default)
  plot_labeled <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", label = TRUE)
  expect_true("GeomText" %in% sapply(plot_labeled$layers, function(l) class(l$geom)[1]))
  
  # Test with label = FALSE
  plot_unlabeled <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", label = FALSE)
  expect_false("GeomText" %in% sapply(plot_unlabeled$layers, function(l) class(l$geom)[1]))
})

test_that("`na.include` parameter works correctly", {
  # Test with na.include = FALSE (default)
  table_no_na <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", exportTable = TRUE, na.include = FALSE)
  expect_false(any(is.na(table_no_na$cloneSize)))
  
  # Test with na.include = TRUE
  table_with_na <- clonalOccupy(scRep_example, x.axis = "seurat_clusters", exportTable = TRUE, na.include = TRUE)
  expect_true(any(is.na(table_with_na$cloneSize)))
  na_row <- table_with_na[table_with_na$seurat_clusters == "1" & is.na(table_with_na$cloneSize), ]
  expect_equal(nrow(na_row), 1)
  expect_equal(na_row$n, 7)
})
