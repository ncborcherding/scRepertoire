# test script for clonalBias.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L"))
scRep_example <- suppressMessages(combineExpression(combined, scRep_example, cloneCall = "aa"))
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Output formats (table vs. plot) are correct", {
  table_output <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                             group.by = "seurat_clusters", n.boots = 2, min.expand = 2, 
                             exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  plot_output <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                            group.by = "seurat_clusters", n.boots = 2, min.expand = 2)
  expect_s3_class(plot_output, "ggplot")
})

test_that("Data frame output has correct structure", {
  table_output <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                             group.by = "seurat_clusters", n.boots = 2, min.expand = 2, 
                             exportTable = TRUE)
  
  expected_cols <- c("Sample", "Clone", "ncells", "Top_state", "bias", "Z.score", "cloneSize")
  expect_true(all(expected_cols %in% names(table_output)))
  expect_true(is.numeric(table_output$ncells))
  expect_true(is.numeric(table_output$bias))
  expect_true(is.numeric(table_output$Z.score))
})

test_that("`min.expand` parameter filters correctly", {
  table_min2 <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                           group.by = "seurat_clusters", n.boots = 2, min.expand = 2, 
                           exportTable = TRUE)
  expect_true(all(table_min2$ncells >= 2))
  table_min1 <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                           group.by = "seurat_clusters", n.boots = 2, min.expand = 1, 
                           exportTable = TRUE)
  expect_gt(nrow(table_min1), nrow(table_min2))
})

test_that("Plot structure is correct", {
  plot_obj <- clonalBias(scRep_example, cloneCall = "aa", split.by = "Patient", 
                         group.by = "seurat_clusters", n.boots = 2, min.expand = 2)
  expect_s3_class(plot_obj$layers[[1]]$geom, "GeomPoint")
  expect_s3_class(plot_obj$layers[[2]]$stat, "StatQuantile")
  expect_equal(plot_obj$labels$x, "Clonal Size")
  expect_equal(plot_obj$labels$y, "Clonal Bias")
  expect_equal(plot_obj$labels$fill, "Group")
  expect_equal(plot_obj$labels$size, "cloneSize")
})

