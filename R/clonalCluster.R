#' Clustering adaptive receptor sequences by edit distance
#'
#' This function uses edit distances of either the nucleotide or amino acid 
#' sequences of the CDR3 and V genes to cluster similar TCR/BCRs together. 
#' As a default, the function takes the input from [combineTCR()], 
#' [combineBCR()] or [combineExpression()] and amends a 
#' cluster to the data frame or meta data. If **exportGraph** is set 
#' to TRUE, the function returns an igraph object of the connected sequences. 
#' If multiple sequences per chain are present, this function only compares
#' the first sequence.
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' sub_combined <- clonalCluster(combined[c(1,2)], 
#'                               chain = "TRA", 
#'                               sequence = "aa")
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()] or [combineExpression()].
#' @param chain Indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param sequence Clustering based on either **"aa"** or 
#' **"nt"**.
#' @param samples The specific samples to isolate for visualization.
#' @param threshold The normalized edit distance to consider. 
#' The higher the number the more similarity of sequence will be 
#' used for clustering.
#' @param group.by The column header used for to group contigs.
#' If (**NULL**), clusters will be calculated across samples.
#' @param use.v,use.j Logical; require identical V/J when `TRUE`
#' @param exportGraph Return an igraph object of connected 
#' sequences (**TRUE**) or the amended input with a
#' new cluster-based variable (**FALSE**).
#' @importFrom igraph set_vertex_attr V union
#' @importFrom plyr join
#' @importFrom stringr str_replace_all
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom stats na.omit
#' @importFrom S4Vectors DataFrame
#' 
#' @export
#' @concept Visualizing_Clones
#' @return Either amended input with edit-distanced clusters added 
#' or igraph object of connect sequences

clonalCluster <- function(input.data, 
                          chain = "TRB", 
                          sequence = "aa",
                          samples = NULL,
                          threshold = 0.85, 
                          group.by = NULL, 
                          use.V = TRUE,
                          use.J = FALSE,
                          exportAdjMatrix = FALSE,
                          exportGraph = FALSE) {
  
  #Prepping any single-cell object
  is.list <- inherits(input.data, "list")
  if(!is.list && !.is.seurat.or.se.object(input.data)) {
    input.data <- .checkList(input.data)
  }
  #TODO Handle NAs
  if(chain == "both") {
    chains <- c("TRA", "TRB")
    test <- lapply(chains, function(x) {
      getIR(input.data, chains = x, sequence.type = sequence, group.by = group.by)
    })
    if (is.data.frame(test[[1]])) {
      detected_chain <- lapply(test, function(df) {
        names(sort(table(substring(df$v, 1, 3))))[1]
      })
    } else {
      detected_chain <- lapply(test, function(inner_list) {
        lapply(inner_list, function(df) {
          names(sort(table(substring(df$v, 1, 3))))[1]
        })
      })
    }
  } else {
    test <- getIR(input.data, chains = chain, sequence.type = sequence)
  }
  
  if (is.data.frame(test[[1]])) {
    network_list <- lapply(test, function(df) {
          tmp <- buildNetwork(df,
                              seq_col   = "cdr3_aa",
                              v_col     = "v",
                              j_col     = "j",
                              filter.v  = use.V,
                              filter.j  = use.J,
                              threshold = threshold)
       
          tmp$from <- df$barcode[as.numeric(tmp$from)]
          tmp$to <- df$barcode[as.numeric(tmp$to)]
          if (threshold < 1) {
            tmp <- tmp[tmp$dist >= threshold, ]
          } else if (threshold %% 1 == 0) {
            tmp <- tmp[tmp$dist <= threshold, ]
          }
        tmp
    })
  } else {
    network_list <- lapply(test, function(inner_list) {
      larger.tmp <- lapply(inner_list, function(df) {
                  tmp <- buildNetwork(df,
                                      seq_col   = "cdr3_aa",
                                      v_col     = "v",
                                      j_col     = "j",
                                      filter.v  = use.V,
                                      filter.j  = use.J,
                                      threshold = threshold)
                  tmp$from <- df$barcode[match(tmp$from, rownames(df))]
                  tmp$to   <- df$barcode[match(tmp$to, rownames(df))]
                  if (threshold < 1) {
                    tmp <- tmp[tmp$dist >= threshold, ]
                  } else if (threshold %% 1 == 0) {
                    tmp <- tmp[tmp$dist <= threshold, ]
                  }
                  tmp
      })
    })
  }
  
  if (is.list(network_list[[1]]) && !is.data.frame(network_list[[1]])) {
    flat_list <- unlist(network_list, recursive = FALSE)
    final_df <- do.call(rbind, flat_list)
  } else {
    final_df <- do.call(rbind, network_list)
  }
  
  graph <- igraph::graph_from_edgelist(as.matrix(final_df[,1:2]))
  
  #Returning the igraph object if exportGraph = TRUE
  if(exportGraph) {
    output.list <- output.list[lapply(output.list,length)>0]
    cluster <- do.call(igraph::union, output.list)
    if(length(is.null(cluster)) == length(cluster)) {
      stop("No clusters detected with current parameters.")
    } else {
      vertex <- names(V(cluster))
      data_df <- unique(data.frame(
        id = vertex
      ))
      data_df <- merge(data_df, graph.variables, by = 1)
      cluster <- set_vertex_attr(cluster, 
                                 name = "size", 
                                 index = data_df$id, 
                                 value = data_df[,2])
      if(ncol(data_df) == 3) { #add grouping variable
        cluster <- set_vertex_attr(cluster, 
                                   name = "group", 
                                   index = data_df$id, 
                                   value = data_df[,3])
      }
      return(cluster)
    }
  }
  
  cluster.list <- lapply(seq_len(length(output.list)), function(x) {
                        output.list[[x]][,1] <- str_replace_all(output.list[[x]][,1], "CTgene", chain)
                        if(!is.null(group.by)) {
                          output.list[[x]][,1] <- paste0(names(dat)[x], ":", output.list[[x]][,1])
                        }
                        colnames(output.list[[x]]) <- c(paste0(chain, "_cluster"), ref2)
                        output.list[[x]]

  })
  cluster <- bind_rows(cluster.list) # the TRA_cluster isnt assigned in the failing test
  
  #Merging with contig info
  tmp <- bound
  tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = TRUE)[,1]
  output2 <- cluster[cluster[,2] %in% tmp[,ref2],]
  bound <-  suppressMessages(join(tmp,  output2, match = "first"))
  
  #Removing unconnected calls
  bound[!grepl(".Cluster", bound[,colnames(cluster)[1]]), colnames(cluster)[1]] <- NA
  
  #Adding to potential single-cell object
  if(inherits(x=input.data, what ="Seurat") | inherits(x=input.data, what ="SummarizedExperiment")) {
    PreMeta <- bound
    x <- colnames(PreMeta)[ncol(PreMeta)]
    PreMeta <- as.data.frame(PreMeta[,x], row.names = rownames(bound))
    colnames(PreMeta) <- x
    if (inherits(x=input.data, what ="Seurat")) { 
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
    #Reorder columns
    bound <- bound[,c(column.order, colnames(cluster)[1])]
    #Split back into list
    if(is.list) {
      if(!is.null(group.by)) {
        bound <- split(bound, bound[,group.by])
      } else {
        bound <- split(bound, bound[,"sample"])
      }
    }
    input.data <- bound
  }
  return(input.data)
}  


  
