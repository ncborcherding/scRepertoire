#' Clustering adaptive receptor sequences
#'
#' This function uses edit distances of either the nucleotide or amino acid 
#' sequences of the CDR3 and V genes to cluster similar TCR/BCRs together. 
#' As a default, the function takes the input from \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}} or \code{\link{combineExpression}} and amends a 
#' cluster to the data frame or meta data. If **exportGraph** is set to TRUE, 
#' the function returns an igraph object of the connected sequences. 
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' sub_combined <- clonalCluster(combined[[2]], chain = "TRA", sequence = "aa")
#' 
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}} 
#' or \code{\link{combineExpression}}.
#' @param chain Indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param sequence Clustering based on either "aa" or "nt".
#' @param samples The specific samples to isolate for visualization.
#' @param threshold The normalized edit distance to consider. The higher the number the more 
#' similarity of sequence will be used for clustering.
#' @param group.by The column header used for to group contigs.
#' @param exportGraph Return an igraph object of connected sequences (TRUE) or the amended
#' input with a new cluster-based variable (FALSE)
#' @importFrom stringdist stringdist
#' @importFrom igraph set_vertex_attr V
#' @importFrom plyr join
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split str_replace_all
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom stats na.omit
#' @export
#' @concept Visualizing_Clones
#' @return Either amended input with edit-distanced clusters added 
#' or igraph object of connect sequences

clonalCluster <- function(df, 
                          chain = "TRB", 
                          sequence = "aa",
                          samples = NULL,
                          threshold = 0.85, 
                          group.by = NULL, 
                          exportGraph = FALSE) {
  output.list <- list()
  if(chain == "both") {
    stop("Please select an individual chain for clustering")
  }
  #Prepping any single-cell object
  is.list <- inherits(df, "list")
  if (inherits(x = df, what = c("Seurat", "SummarizedExperiment"))) {
    dat <- .grabMeta(df)
  } else {
    dat <- df
  }
  dat <- .checkList(dat)
  dat <- .data.wrangle(dat, group.by, "CTgene", chain)
  
  if (inherits(x = df, what = c("Seurat", "SummarizedExperiment"))) {
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
  dat <- .checkList(dat)
  
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
    
    graph.variables <- bound %>%
                          group_by(bound[,ref2]) %>%
                          summarize(sample_count = n(),
                                    unique_samples = paste0(unique(group.by), collapse = ","))
    
  } else {
    bound <- bind_rows(dat)
    graph.variables <- bound %>%
        group_by(bound[,ref2]) %>%
        summarize(sample_count = n())
  }

  #Generating Connected Component
  cluster <- .lvCompare(bound, 
                        gene = "CTgene", 
                        chain = ref2, 
                        threshold = threshold, 
                        exportGraph = exportGraph)
  
  if(exportGraph) {
    vertex <- names(V(cluster))
    data_df <- unique(data.frame(
      id = V(cluster)$name
    ))
    data_df <- merge(data_df, graph.variables, by = 1)
    cluster <- set_vertex_attr(cluster, name = "size", index = data_df$id, value = data_df[,2])
    if(ncol(data_df) == 3) { #add grouping variable
      cluster <- set_vertex_attr(cluster, name = "group", index = data_df$id, value = data_df[,3])
    }
    return(cluster)
  }
  
  cluster[,1] <-  str_replace_all(cluster[,1], "CTgene", chain)
  colnames(cluster) <- c(paste0(chain, "_cluster"), ref2)
  
  #Merging with contig info
  tmp <- bound
  tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = TRUE)[,1]
  output2 <- cluster[cluster[,2] %in% tmp[,ref2],]
  bound <-  suppressMessages(join(tmp,  output2, match = "first"))
  
  #Adding to potential single-cell object
  if(inherits(x=df, what ="Seurat") | inherits(x=df, what ="SummarizedExperiment")) {
    PreMeta <- bound
    x <- colnames(PreMeta)[ncol(PreMeta)]
    PreMeta <- as.data.frame(PreMeta[,x], row.names = rownames(bound))
    colnames(PreMeta) <- x
    if (inherits(x=df, what ="Seurat")) { 
      col.name <- names(PreMeta) %||% colnames(PreMeta)
      df[[col.name]] <- PreMeta
    } else {
      rownames <- rownames(colData(df))
      colData(df) <- cbind(colData(df), PreMeta[rownames,])[, union(colnames(colData(df)),  colnames(PreMeta))]
      rownames(colData(df)) <- rownames 
    }
  } else {
    df <- bound
  }
  return(df)
}  


  