# test script for clonalCluster.R - testcases are NOT comprehensive!

# Data to Use
combined <- combineTCR(contig_list,
                      samples = c("P17B", "P17L", "P18B", "P18L",
                                  "P19B", "P19L", "P20B", "P20L"))


test_that("Basic functionality and default output structure", {
  clustered_list <- clonalCluster(combined)
  expect_length(clustered_list, length(combined))
  expect_true(all(sapply(clustered_list, is.data.frame)))
  expect_true("TRB.Cluster" %in% names(clustered_list[[1]]))
  
  # Check for data integrity
  expect_equal(nrow(clustered_list[[1]]), nrow(combined[[1]]))
  expect_false(all(is.na(clustered_list[[1]]$TRB.Cluster)))
})


test_that("exportGraph = TRUE returns a valid igraph object", {
  graph_obj <- clonalCluster(combined, exportGraph = TRUE)
  expect_s3_class(graph_obj, "igraph")
  expect_gt(vcount(graph_obj), 0)
  expect_gt(ecount(graph_obj), 0)
  expect_true("cluster" %in% igraph::vertex_attr_names(graph_obj))
  expect_true("weight" %in% igraph::edge_attr_names(graph_obj))
})


test_that("exportAdjMatrix = TRUE returns a valid sparse matrix", {
  adj_matrix <- clonalCluster(combined, exportAdjMatrix = TRUE)
  all_barcodes <- unique(do.call(rbind, combined)[["barcode"]])
  num_barcodes <- length(all_barcodes)
  expect_s4_class(adj_matrix, "dgCMatrix")
  expect_equal(dim(adj_matrix), c(num_barcodes, num_barcodes))
  expect_equal(rownames(adj_matrix), all_barcodes)
})


test_that("chain parameter works correctly", {
  clustered_tra <- clonalCluster(combined, chain = "TRA")
  expect_true("TRA.Cluster" %in% names(clustered_tra[[1]]))
  expect_false("TRB.Cluster" %in% names(clustered_tra[[1]]))
  clustered_both <- clonalCluster(combined, chain = "both")
  expect_true("Multi.Cluster" %in% names(clustered_both[[1]]))
})


test_that("group.by parameter functions without error", {
  # Grouping by sample should produce clusters within each sample
  clustered_grouped <- clonalCluster(combined, group.by = "Sample")
  
  # Basic structural checks
  expect_s3_class(clustered_grouped, "list")
  expect_true("TRB_Cluster" %in% names(clustered_grouped[[1]]))
  
  # A simple check to ensure clustering happened
  # More complex checks would require knowing the expected ground truth
  expect_false(all(is.na(clustered_grouped[[1]]$TRB_Cluster)))
})


test_that("threshold parameter influences clustering", {
  lenient_clusters <- clonalCluster(combined, threshold = 5)
  num_lenient <- sum(!is.na(lenient_clusters[[7]]$TRB.Cluster))
  strict_clusters <- clonalCluster(combined, threshold = 0.95)
  num_strict <- sum(!is.na(strict_clusters[[7]]$TRB.Cluster))
  
  # Expect more vertices to be part of a cluster with a lenient threshold
  expect_gt(num_lenient, num_strict)
})


test_that("Different `cluster.method` options work", {
  walktrap_graph <- clonalCluster(combined, cluster.method = "walktrap", exportGraph = TRUE)
  expect_s3_class(walktrap_graph, "igraph")
  expect_true("cluster" %in% igraph::vertex_attr_names(walktrap_graph))
})


test_that("8. Input validation and error handling", {
  # Error when both export options are TRUE
  expect_error(
    clonalCluster(combined, exportGraph = TRUE, exportAdjMatrix = TRUE),
    "Please set only one of `exportGraph` or `exportAdjMatrix` to TRUE."
  )
  
  # Error on invalid cluster method
  expect_error(
    clonalCluster(combined, cluster.method = "invalid_method"),
    "Unsupported clustering.method: 'invalid_method'"
  )
})



