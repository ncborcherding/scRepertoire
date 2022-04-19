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
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' sub_combined <- clusterTCR(combined[[2]], chain = "TRA", sequence = "aa")
#' 
#' @param df The product of combineTCR(), expression2List(), or combineExpression().
#' @param chain The TCR to cluster - TRA, TRB, TRG, TRD
#' @param sequence Clustering based on either "aa" or "nt"
#' @param threshold The normalized edit distance to consider. The higher the number the more 
#' similarity of sequence will be used for clustering.
#' @param group.by The column header used for to calculate the cluster
#' @importFrom stringdist stringdistmatrix
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
    
    for (x in seq_along(bound)) {
      dictionary <- na.omit(unique(bound[[x]][,ref2]))
      dictionary <- str_split(dictionary, ";", simplify = TRUE)[,1]
      length <- nchar(dictionary)
      dist <- stringdistmatrix(dictionary, method = "lv")  
      edge.list <- NULL
      for (j in seq_len(length(dictionary))) {
        row <- SliceExtract_dist(dist,j)
        norm.row <- row
        for (k in seq_len(length(norm.row))) {
          norm.row[k] <- 1- (norm.row[k]/mean(c(length[j],length[k])))
        }
        neighbor <- which(norm.row >= threshold)
        knn.norm = data.frame("from" = j,
                              "to" = neighbor)
        edge.list <- rbind(edge.list, knn.norm)
      }
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
      output.list[[x]] <- output
    }
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

#Code from https://stackoverflow.com/questions/57282842/how-to-efficiently-extract-a-row-or-column-from-a-dist-distance-matrix?rq=1
f <- function (i, j, dist_obj) {
      if (!inherits(dist_obj, "dist")) stop("please provide a 'dist' object")
      n <- attr(dist_obj, "Size")
      valid <- (i >= 1) & (j >= 1) & (i > j) & (i <= n) & (j <= n)
      k <- (2 * n - j) * (j - 1) / 2 + (i - j)
      k[!valid] <- NA_real_
      k
}
    
#Code from https://stackoverflow.com/questions/57282842/how-to-efficiently-extract-a-row-or-column-from-a-dist-distance-matrix?rq=1
SliceExtract_dist <- function (dist_obj, k) {
      if (length(k) > 1) stop("The function is not 'vectorized'!")
      n <- attr(dist_obj, "Size")
      if (k < 1 || k > n) stop("k out of bound!")
      ##
      i <- 1:(k - 1)
      j <- rep.int(k, k - 1)
      v1 <- dist_obj[f(j, i, dist_obj)]
      ## 
      i <- (k + 1):n
      j <- rep.int(k, n - k)
      v2 <- dist_obj[f(i, j, dist_obj)]
      ## 
      c(v1, 0, v2)
}
