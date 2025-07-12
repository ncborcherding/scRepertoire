# test script for clonalNetwork.R - testcases are NOT comprehensive!
library(ggraph)

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)
 
test_that("Output formats (table, clones, plot) are correct", {
  # Test exportTable
  table_output <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  expect_true(all(c("to", "from", "weight") %in% names(table_output)))
  
  # Test exportClones
  clones_output <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", exportClones = TRUE)
  expect_s3_class(clones_output, "data.frame")
  expect_true(all(c("clone", "sum") %in% names(clones_output)))
  
  # Test plot output
  plot_output <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters")
  expect_s3_class(plot_output, "ggplot")
})

test_that("Core edge weight calculation is correct", {
  edge_list <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", exportTable = TRUE)
  
  edge_c3_c5 <- edge_list[edge_list$from == "3" & edge_list$to == "5", ]
  edge_c5_c3 <- edge_list[edge_list$from == "5" & edge_list$to == "3", ]
  
  expect_equal(edge_c3_c5$weight, 0.118, tolerance = 1e-2)
  expect_equal(edge_c5_c3$weight, 0.1, tolerance = 1e-2)
})

test_that("Filtering parameters work as expected", {
  edge_list_prop <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", 
                                  filter.proportion = 0.2, exportTable = TRUE)
  expect_equal(nrow(edge_list_prop), 1)
  expect_equal(edge_list_prop$weight, 0.25)
  edge_list_unfiltered <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", exportTable = TRUE)
  edge_list_half <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", 
                                  filter.graph = TRUE, exportTable = TRUE)
  expect_lt(nrow(edge_list_half), nrow(edge_list_unfiltered))
})

test_that("exportClones provides correctly summarized data", {
  clones_df <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters", exportClones = TRUE)
  expect_equal(nrow(clones_df), 341)
  expect_equal(clones_df$sum[1:5], c(11,3,3,3,2))
})

test_that("Plot structure has correct geoms", {
  plot_obj <- clonalNetwork(scRep_example, reduction = "umap", group.by = "seurat_clusters")
  geoms_in_plot <- sapply(plot_obj$layers, function(l) class(l$geom)[1])
  expect_true("GeomPoint" %in% geoms_in_plot)    
  expect_true("GeomEdgePath" %in% geoms_in_plot)  
})


