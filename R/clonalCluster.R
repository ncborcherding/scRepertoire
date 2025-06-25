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
#' @param exportGraph Return an igraph object of connected 
#' sequences (**TRUE**) or the amended input with a
#' new cluster-based variable (**FALSE**).
#' @importFrom stringdist stringdist
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
                          exportGraph = FALSE) {
  if(chain == "both") {
    stop("Please select an individual chain for clustering")
  }
  #Prepping any single-cell object
  is.list <- inherits(input.data, "list")
  if(!is.list && !inherits(x = input.data, what = c("Seurat", "SummarizedExperiment"))) {
    input.data <- .checkList(input.data)
  }
  dat <- .dataWrangle(input.data, group.by, "CTgene", chain)
  
  if (inherits(x = input.data, what = c("Seurat", "SummarizedExperiment"))) {
    for (y in seq_along(dat)) {
      cdr3_aa <- str_split(dat[[y]]$CTaa, "_", simplify = TRUE)
      cdr3_nt <- str_split(dat[[y]]$CTnt, "_", simplify = TRUE)
      
      dat[[y]]$cdr3_aa1 <- cdr3_aa[, 1]
      dat[[y]]$cdr3_aa2 <- cdr3_aa[, 2]
      dat[[y]]$cdr3_nt1 <- cdr3_nt[, 1]
      dat[[y]]$cdr3_nt2 <- cdr3_nt[, 2]
    }
  }
  if (!is.null(samples)) {
    dat <- dat[which(names(dat) %in% samples)]
  }
  
  #Getting column labels to pull from
  if (grepl("TRA|TRG", chain)) {
    gene <- "TCR"
    ref <- 1
  } else if (grepl("TRB|TRD", chain)) {
    gene <- "TCR"
    ref <- 2
  } else if (grepl("IGH", chain)) {
    gene <- "BCR"
    ref <- 1
  } else if (grepl("IGL", chain)) {
    gene <- "BCR"
    ref <- 2
  }
  ref2 <- paste0("cdr3_", sequence, ref)
  
  if (!is.null(group.by)) {
    bound <- bind_rows(dat, .id = "group.by")
    bound <- bound[!is.na(bound[,ref2]),]
    retain.ref <- data.frame(old = bound[,ref2], new = str_split(bound[,ref2], ";", simplify = TRUE)[,1])
    bound[,ref2] <- str_split(bound[,ref2], ";", simplify = TRUE)[,1]
    graph.variables <- bound %>%
                          group_by(bound[,ref2]) %>%
                          dplyr::summarize(sample_count = n(),
                                    unique_samples = paste0(unique(group.by), collapse = ","))
  } else {
    bound <- bind_rows(dat)
    bound <- bound[!is.na(bound[,ref2]),]
    retain.ref <- data.frame(old = bound[,ref2], new = str_split(bound[,ref2], ";", simplify = TRUE)[,1])
    bound[,ref2] <- str_split(bound[,ref2], ";", simplify = TRUE)[,1]
    graph.variables <- bound %>%
                          group_by(bound[,ref2]) %>%
                          dplyr::summarize(sample_count = n())
  }
  dictionary <- dat
  #Generating Connected Component
  output.list <- lapply(dictionary, function(x) {
    cluster <- .lvCompare(x, 
                          gene = "CTgene", 
                          chain = ref2, 
                          threshold = threshold,  
                          exportGraph = exportGraph)
    cluster
  })
  #Grabbing column order for later return
  column.order <- colnames(bound)
  
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


  
