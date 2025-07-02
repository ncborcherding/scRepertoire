#' Cluster Clones by sequence similarity
#'
#' This function clusters TCRs or BCRs based on the edit distance of their CDR3
#' sequences. It can operate on either nucleotide ("nt") or amino acid ("aa")
#' sequences and can optionally enforce that clones share the same V gene.
#' The output can be the input object with an added metadata column for cluster
#' IDs, a sparse adjacency matrix, or an `igraph` graph object representing the
#' cluster network.
#' 
#' @details
#' The clustering process is as follows:
#' 1.  The function retrieves the relevant chain data from the input object.
#' 2.  It calculates the edit distance between all sequences within each group
#'     (or across the entire dataset if `group.by` is `NULL`).
#' 3.  An edge list is constructed, connecting sequences that meet the similarity
#'     `threshold`.
#' 4.  The `threshold` parameter behaves differently based on its value:
#'     - **`threshold` < 1 (e.g., 0.85):** Interpreted as a *normalized* edit
#'       distance or sequence similarity. A higher value means greater
#'       similarity is required. This is the default behavior.
#'     - **`threshold` >= 1 (e.g., 2):** Interpreted as a maximum *raw* edit
#'       distance. A lower value means greater similarity is required.
#' 5.  An `igraph` graph is built from the edge list.
#' 6.  A clustering algorithm is run on the graph. The default 
#'     `cluster.method = "components"` simply identifies the connected 
#'     components (i.e., each cluster is a group of sequences connected by 
#'     edges). Other methods from `igraph` can be used.
#' 7.  The resulting cluster information is formatted and returned in the
#'     specified format.
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B","P19L", "P20B", "P20L"))
#'
#' # Add cluster information to the list
#' sub_combined <- clonalCluster(combined[c(1,2)],
#'                               chain = "TRA",
#'                               sequence = "aa",
#'                               threshold = 0.85)
#'
#' # Export the graph object instead
#' graph_obj <- clonalCluster(combined[c(1,2)],
#'                            chain = "TRA",
#'                            exportGraph = TRUE)
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()] or [combineExpression()].
#' @param chain The TCR/BCR chain to use for clustering. Use "both" to combine
#' networks from TRA/TRB or IGH/IGL. Valid options: "both", "TRA", "TRB",
#' "TRD", "TRG", "IGH", "IGK", "IGL".
#' @param sequence Clustering based on either **"aa"** or 
#' **"nt"**.
#' @param threshold The similarity threshold. If < 1, treated as normalized
#' similarity (higher is stricter). If >= 1, treated as raw edit distance
#' (lower is stricter).
#' @param group.by The column header used for to group clones.
#' If (**NULL**), clusters will be calculated across samples.
#' @param use.V Logical. If `TRUE`, sequences must share the same V gene to be
#' clustered together.
#' @param use.J Logical. If `TRUE`, sequences must share the same J gene to be
#' clustered together.
#' @param cluster.method The clustering algorithm to use. Defaults to `"components"`, 
#' which finds connected subgraphs.
#' @param cluster.prefix A character prefix to add to the cluster names (e.g.,
#' "cluster.").
#' @param exportGraph Logical. If `TRUE`, the function returns an `igraph`
#' object of the sequence network.
#' @param exportAdjMatrix Logical. If `TRUE`, the function returns a sparse
#' adjacency matrix (`dgCMatrix`) of the network.
#' @param exportGraph Return an igraph object of connected 
#' sequences (**TRUE**) or the amended input with a
#' new cluster-based variable (**FALSE**).
#' @importFrom igraph graph_from_edgelist E E<- V V<- as_data_frame 
#' as_adjacency_matrix membership
#' @importFrom dplyr left_join
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom S4Vectors DataFrame
#' @importFrom Matrix sparseMatrix
#' 
#' @export
#' @concept Visualizing_Clones
#' @return 
#' Depending on the export parameters, one of the following:
#' \itemize{
#'   \item An amended `input.data` object with a new metadata column containing cluster IDs (default).
#'   \item An `igraph` object if `exportGraph = TRUE`.
#'   \item A sparse `dgCMatrix` object if `exportAdjMatrix = TRUE`.
#' }

clonalCluster <- function(input.data, 
                          chain = "TRB", 
                          sequence = "aa",
                          samples = NULL,
                          threshold = 0.85, 
                          group.by = NULL, 
                          cluster.method = "components",
                          cluster.prefix = "cluster.",
                          use.V = TRUE,
                          use.J = FALSE,
                          exportAdjMatrix = FALSE,
                          exportGraph = FALSE) {
  
  if (exportGraph && exportAdjMatrix) {
    stop("Please set only one of `exportGraph` or `exportAdjMatrix` to TRUE.")
  }
  
  #Prepping any single-cell object
  is.list <- inherits(input.data, "list")
  if(!is.list && !.is.seurat.or.se.object(input.data)) {
    input.data <- .checkList(input.data)
  }
  if (chain == "both") {
    # Use a placeholder to grab both alpha/beta or heavy/light chains
    chains_to_get <- c("TRA", "TRB")
    chain_data <- lapply(chains_to_get, function(x) {
      getIR(input.data, chains = x, sequence.type = sequence, group.by = group.by)
    })
  } else {
    chain_data <- getIR(input.data, chains = chain, sequence.type = sequence, group.by = group.by)
    chain_data <- list(chain_data)
  }
  
  # Flatten the list if nested and get all barcodes
  if (is.list(chain_data[[1]]) && !is.data.frame(chain_data[[1]])) {
    chain_data <- unlist(chain_data, recursive = FALSE)
  }
  all_barcodes <- unique(do.call(rbind, chain_data)[["barcode"]])
  
  # Apply the network function to each data frame and combine into one edge list
  full_edge_list <- do.call(rbind, lapply(chain_data, .buildNetwork))
  full_edge_list <- unique(full_edge_list) 
  
  if (nrow(full_edge_list) == 0) {
    warning("No clusters found, returning original data.")
    return(input.data)
  }
  
  # Create the graph object
  full_g <- igraph::graph_from_edgelist(as.matrix(full_edge_list[, c("from", "to")]), directed = FALSE)
  igraph::E(full_g)$weight <- full_edge_list$dist
  
  # Perform clustering using the specified method
  clusters <- .clusterGraph(cluster.method, full_g)
  cluster_membership <- igraph::membership(clusters)
  
  # Rename clusters sequentially with the specified prefix
  unique_clusters <- unique(cluster_membership)
  cluster_names <- paste0(cluster.prefix, seq_along(unique_clusters))
  name_map <- setNames(cluster_names, unique_clusters)
  renamed_membership <- name_map[as.character(cluster_membership)]
  names(renamed_membership) <- names(cluster_membership)
  
  # Add vertex communities 
  igraph::V(full_g)$cluster <- renamed_membership
  
  # Returning Graph
  if(exportGraph) {
    # Adding additional vertex information if graph being exported
    vertex_df <- igraph::as_data_frame(full_g, what = "vertices")
    colnames(vertex_df)[1] <- "barcode" 
    full_meta_long <- do.call(rbind, chain_data)
    
    if (nrow(full_meta_long) > 0) {
      meta_to_process <- unique(full_meta_long[, c("barcode", "cdr3_aa", "v", "j")])
      
      # Collapsing incase > 1 chain
      meta_indexed <- meta_to_process %>% 
                          group_by(barcode) %>%
                          mutate(chain_num = dplyr::row_number()) %>%
                          ungroup() %>%
                          as.data.frame()
      original_meta <- reshape(
        data = meta_indexed, 
        idvar = "barcode",
        timevar = "chain_num",
        direction = "wide",
        v.names = c("cdr3_aa", "v", "j"),
        sep = "" )
      # Loop to add more vertex info
      vertex_names <- V(full_g)$name
      match_indices <- match(vertex_names, original_meta$barcode)
      for (col in colnames(original_meta)) {
        if (col == "barcode") next 
        full_g <- set_vertex_attr(
          graph = full_g,
          name = col,
          value = original_meta[[col]][match_indices])
      }
    }
    return(full_g)
  }
  
  # Returning Adjacency Matrix
  if (exportAdjMatrix) {
    adj_from_graph <- igraph::as_adjacency_matrix(full_g, 
                                                  sparse = TRUE,
                                                  attr = "weight")
    barcodes_in_graph <- rownames(adj_from_graph)
    
    adjacency_matrix <- Matrix::sparseMatrix(
      i = {}, j = {}, x = 0.0,
      dims = c(length(all_barcodes), length(all_barcodes)),
      dimnames = list(all_barcodes, all_barcodes)
    )
    adjacency_matrix[barcodes_in_graph, barcodes_in_graph] <- adj_from_graph
    return(adjacency_matrix)
  }

 # Attaching to input.data
 bound <- igraph::as_data_frame(full_g, what = "vertices")
 colnames(bound)[2] <- ifelse(chain == "both", "Multi.Cluster", paste0(chain, ".Cluster"))
  
  #Adding to potential single-cell object
  if(.is.seurat.or.se.object(input.data)) {
    PreMeta <- bound[,-1, drop = FALSE]
    if (.is_seurat_object(input.data)) { 
      col.name <- names(PreMeta) %||% colnames(PreMeta)
      input.data[[col.name]] <- PreMeta
    } else {
      combined_col_names <- unique(c(colnames(colData(sc.data)), colnames(PreMeta)))
      full_data <- merge(colData(sc.data), PreMeta[rownames, , drop = FALSE], by = "row.names", all.x = TRUE)
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(sc.data) <- DataFrame(full_data[, combined_col_names])
    }
  } else {
    colnames(bound)[1] <- "barcode"
    input.data <- lapply(input.data, function(df) {
      dplyr::left_join(df, bound, by = "barcode")
    })
  }
  return(input.data)
}  

# Handling clustering options
#' @importFrom igraph cluster_fast_greedy cluster_walktrap cluster_louvain 
#' cluster_leiden cluster_spinglass cluster_edge_betweenness components
clusterGraph <- function(clustering.method, graph, ...) {
  if (!igraph::is_igraph(graph)) {
    stop("The 'graph' argument must be a valid 'igraph' object.", call. = FALSE)
  }
  
  if (!is.character(clustering.method) || length(clustering.method) != 1) {
    stop("The 'clustering.method' argument must be a single character string.", call. = FALSE)
  }
  result <- switch(
    clustering.method,
    "fast_greedy" = igraph::cluster_fast_greedy(graph, ...),
    "walktrap" = igraph::cluster_walktrap(graph, ...),
    "louvain" = igraph::cluster_louvain(graph, ...),
    "leiden" = igraph::cluster_leiden(graph, ...),
    "spinglass" = igraph::cluster_spinglass(graph, ...),
    "edge_betweenness" = igraph::cluster_edge_betweenness(graph, ...),
    "components" = igraph::components(graph, ...),
    
    # Default case: Handle unsupported methods
    stop(
      paste0(
        "Unsupported clustering.method: '", clustering.method, "'. ",
        "Supported methods are: 'fast_greedy', 'walktrap', 'louvain', ",
        "'leiden', 'spinglass', 'edge_betweenness', and 'components'."
      ),
      call. = FALSE
    )
  )
  
  return(result)
}

#' @importFrom immApex buildNetwork
.buildNetwork <- function(df) {
  edge_list <- buildNetwork(df,
                            seq_col   = "cdr3_aa",
                            v_col     = "v",
                            j_col     = "j",
                            filter.v  = use.V,
                            filter.j  = use.J,
                            ids = df[["barcode"]],
                            threshold = threshold)
  
  return(edge_list)
}
  
