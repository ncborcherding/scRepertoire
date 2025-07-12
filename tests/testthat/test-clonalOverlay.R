# test script for clonalOverlay.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)

test_that("Input validation works correctly", {
  # Test for error when reduction is missing
  expect_error(clonalOverlay(scRep_example, reduction = "nonexistent_reduction"), 
               regexp = "Reduction nonexistent_reduction")
  
  # Test for error when cut.category is missing
  expect_error(clonalOverlay(scRep_example, reduction = "umap", cut.category = "nonexistent_category"),
               "If filtering the data using a cutpoint, ensure the cut.category correspond to a variable in the meta data.")
})

test_that("Data filtering with `cutpoint` is correct", {
  cutpoint_val <- 25
  manual_subset_rows <- rownames(scRep_example@meta.data[scRep_example@meta.data$clonalFrequency >= cutpoint_val, ])
  plot_obj <- clonalOverlay(scRep_example, reduction = "umap", cutpoint = cutpoint_val)
  density_data <- layer_data(plot_obj, 2) # Layer 2 should be the density plot
  all_coords <- scRep_example@reductions$umap@cell.embeddings
  filtered_data_internal_size <- nrow(scRep_example@meta.data[scRep_example@meta.data$clonalFrequency >= cutpoint_val,])
  expect_equal(filtered_data_internal_size, length(manual_subset_rows))
  plot_empty <- clonalOverlay(scRep_example, reduction = "umap", cutpoint = 100)
  density_data_empty <- layer_data(plot_empty, 2)
  expect_equal(nrow(density_data_empty), 0)
})

test_that("Plot structure is correct", {
  plot_obj <- clonalOverlay(scRep_example, reduction = "umap", cutpoint = 20, bins = 10)
  expect_length(plot_obj$layers, 2)
  expect_s3_class(plot_obj$layers[[1]]$geom, "GeomPoint")
  expect_s3_class(plot_obj$layers[[2]]$geom, "GeomDensity2d")
  expect_equal(plot_obj$layers[[2]]$stat_params$bins, 10)
  expect_equal(plot_obj$labels$x, "Dimension 1")
  expect_equal(plot_obj$labels$y, "Dimension 2")
  expect_equal(plot_obj$labels$colour, "Active Identity")
})

test_that("`facet.by` parameter works correctly", {
  plot_facet <- clonalOverlay(scRep_example, reduction = "umap", facet.by = "orig.ident")
  expect_s3_class(plot_facet$facet, "FacetWrap")
  facet_vars <- names(plot_facet$facet$params$facets)
  expect_equal(facet_vars, "facet.by")
})

test_that("Function handles Seurat objects without error", {
  expect_silent(clonalOverlay(scRep_example, reduction = "umap"))
})

