#' Clustering T cell receptors
#'
#' This function uses edit distances of either the nucleotide or amino acid 
#' sequences of the CDR3 to cluster similar TCRs together. The distance clustering
#' will then be amended to the end of the list of combined contigs. The cluster 
#' will appear as CHAIN.num if a unique sequence or CHAIN:LD.num if clustered together.
#' This function will only two chains recovered, multiple chains will automatically 
#' be reduced. This function also underlies the combineBCR() function and therefore 
#' not needed for B cells. This may take some time to calculate the distances and cluster. 
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
#' @param group.by The column header used for to calculate the cluster
#' @importFrom stringdist stringdist
#' @importFrom igraph graph_from_data_frame components
#' @importFrom plyr join
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split
#' @importFrom  rlang %||%
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
  
  output.list <- lapply(bound, function(x) {
    dictionary <- na.omit(unique(x[,ref2]))
    dictionary <- str_split(dictionary, ";", simplify = TRUE)[,1]
    dictionary <- na.omit(unique(dictionary))
    length <- nchar(dictionary)
    edge.list <- lapply(dictionary, function(y) {
      pos <- which(dictionary == y)
      dist <- stringdist(y, dictionary[-pos], method = "lv")
      norm.row <- dist
      norm.row <- 1 - (norm.row/((length[pos] + length[-pos])/2))
      neighbor <- which(norm.row >= threshold)
      if(length(neighbor) > 0) {
        nn.norm = data.frame("from" = pos,
                             "to" = neighbor)
      } else {
        nn.norm <- NULL
      }
      nn.norm
    })
    edge.list = edge.list[-which(sapply(edge.list, is.null))]
    edge.list <- do.call(rbind, edge.list)
    edge.list$to <-ifelse(edge.list$to > edge.list$from, edge.list$to + 1, edge.list$to)
    edge.list <- unique(edge.list)
    g <- graph_from_data_frame(edge.list)
    components <- components(g, mode = c("weak"))
    out <- data.frame("cluster" = components$membership, 
                      "filtered" = names(components$membership))
    filter <- which(table(out$cluster) > 1)
    out <- subset(out, cluster %in% filter)
    if(nrow(out) > 1) {
      if (list.length == 1) {
        out$cluster <- paste0(chain, ":LD", ".", out$cluster)
      } else {
        out$cluster <- paste0(names(bound)[x], ".", chain, ":LD", ".", out$cluster)
      }
      out$filtered <- dictionary[as.numeric(out$filtered)]
    }
    uni_IG <- as.data.frame(unique(dictionary[dictionary %!in% out$filtered]))
    colnames(uni_IG) <- "filtered"
    if (nrow(uni_IG) > 0) {
      if (list.length == 1) {
        uni_IG$cluster <- paste0(chain, ".", seq_len(nrow(uni_IG))) 
      } else {
        uni_IG$cluster <- paste0(names(bound)[x], ".", chain, ".", seq_len(nrow(uni_IG))) 
      }
    }
    
    output <- rbind.data.frame(out, uni_IG)
    colname <- paste0(chain, "_cluster")
    colnames(output) <- c(colname,ref2)
    output
  })
  for (i in seq_along(bound)) {
    tmp <- bound[[i]]
    output <- bind_rows(output.list)
    tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = TRUE)[,1]
    output2 <- output[output[,2] %in% tmp[,ref2],]
    
    tmp <-  unique(suppressMessages(join(tmp,  output2)))
    bound[[i]] <- tmp
  }
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