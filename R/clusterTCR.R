#' Clustering T cell receptors
#'
#' This function uses edit distances of either the nucleotide or amino acid 
#' sequences of the CDR3 and V genes to cluster similar TCRs together. 
#' The distance clustering will then be amended to the end of the list of 
#' combined contigs. The cluster will appear as CHAIN.num if a unique sequence or 
#' CHAIN:LD.num if clustered together.T his function will only two chains 
#' recovered, multiple chains will automatically be reduced. This function 
#' also underlies the combineBCR() function and therefore not needed for 
#' B cells. This may take some time to calculate the distances and cluster. 
#' 
#' @examples
# Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' 
#' sub_combined <- clusterTCR(combined[[2]], chain = "TRA", sequence = "aa")
#' 
#' @param df The product of combineTCR(), expression2List(), or combineExpression().
#' @param chain The TCR to cluster - TRA, TRB, TRG, TRD
#' @param sequence Clustering based on either "aa" or "nt"
#' @param threshold The normalized edit distance to consider. The higher the number the more 
#' similarity of sequence will be used for clustering.
#' @param group.by The column header used for to group contigs.
#' @importFrom stringdist stringdist
#' @importFrom igraph graph_from_data_frame components
#' @importFrom plyr join
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split str_replace_all
#' @importFrom rlang %||%
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom stats na.omit
#' @export
#' @return List of clonotypes for individual cell barcodes

clusterTCR <- function(df, 
                       chain = NULL, 
                       sequence = NULL, 
                       threshold = 0.85, 
                       group.by = NULL) {
  output.list <- list()
  dat <- list.input.return(df, group.by)
  if (inherits(x=df, what ="Seurat") |
      inherits(x=df, what ="SummarizedExperiment")) {
    for(y in seq_along(dat)) {
      dat[[y]]$cdr3_aa1 <- str_split(dat[[y]]$CTaa, "_", simplify = TRUE)[,1]
      dat[[y]]$cdr3_aa2 <- str_split(dat[[y]]$CTaa, "_", simplify = TRUE)[,2]
      dat[[y]]$cdr3_nt1 <- str_split(dat[[y]]$CTnt, "_", simplify = TRUE)[,1]
      dat[[y]]$cdr3_nt2 <- str_split(dat[[y]]$CTnt, "_", simplify = TRUE)[,2]
    }
  }
  dat <- checkList(dat)
  if(chain %in% c("TRA", "TRG")) {
    ref <- 1
  } else if(chain %in% c("TRB", "TRD")) {
    ref <- 2
  }
  ref2 <- paste0("cdr3_", sequence, ref)
  bound <- bind_rows(dat)
  #Should make it work as either grouped or non-grouped
  if (!is.null(group.by)) {
    bound <- split(bound, bound[,group.by])
    list.length <- length(bound)
  } else {
    bound <- list(bound)
    list.length <- 1
  }
  #Defining clusters by edit distance
  output.list <- lapply(bound, function(x) {
    cluster <- lvCompare(x, gene = paste0("TCR", ref), chain = ref2, threshold = threshold)
    cluster
  })
  #Reformating output for compatibility
  for (i in seq_along(output.list)) {
        output.list[[i]][,1] <- str_replace_all(output.list[[i]][,1], paste0("TCR", ref), chain)
        if (length(output.list) > 1) {
          output.list[[i]][,1] <- paste0(names(bound)[i], ".", output.list[[i]][,1]) 
        }
        colnames(output.list[[i]])[1] <- paste0(chain, "_cluster")
        colnames(output.list[[i]])[2] <- paste0(ref2)
  }
  #Merging with contig info
  for (i in seq_along(bound)) {
    tmp <- bound[[i]]
    output <- bind_rows(output.list)
    tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = TRUE)[,1]
    output2 <- output[output[,2] %in% tmp[,ref2],]
    
    tmp <-  unique(suppressMessages(join(tmp,  output2)))
    bound[[i]] <- tmp
  }
  #Adding to potential single-cell object
  if(inherits(x=df, what ="Seurat") | inherits(x=df, what ="SummarizedExperiment")) {
    PreMeta <- bind_rows(bound)
    x <- colnames(PreMeta)[ncol(PreMeta)]
    PreMeta <- as.data.frame(PreMeta[,x], row.names = PreMeta$barcode)
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

