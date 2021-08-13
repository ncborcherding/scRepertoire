#' Clustering T cell receptors
#'
#' This function uses edit distances of either the nucleotide or amino acid 
#' sequences of the CDR3 to cluster similar TCRs together. The distance clustering
#' will then be ammended to the end of the list of combined contigs with the corresponding Vgene. The
#' cluster will appear as CHAIN.num if a unique sequence or CHAIN:LD.num if clustered together.
#' This function will only two chains recovered, multiple chains will automatically be reduced. This 
#' function also underlies the combineBCR() function and therefore not needed for B cells. This may
#' take some time to calculate the distances and cluster. 
#' 
#' @examples
# Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' sub_combined <- clusterTCR(combined[[2]], chain = "TCRA", sequence = "aa")
#' 
#' @param df The product of CombineTCR() or CombineBCR().
#' @param chain The TCR to cluster
#' @param sequence Clustering based on either "aa" or "nt"
#' @param threshold The normalized edit distance to consider. The higher the number the more 
#' similarity of sequence will be used for clustering.
#' @param group The column header used for to calculate the cluster by.
#' @importFrom stringdist stringdistmatrix
#' @importFrom igraph graph_from_data_frame components
#' @importFrom plyr join
#' @importFrom dplyr bind_rows
#' @export
#' @return List of clonotypes for individual cell barcodes

clusterTCR <- function(df, chain = NULL, sequence = NULL, threshold = 0.85, group = NULL) {
    output.list <- list()
    `%!in%` = Negate(`%in%`)
    df <- checkList(df)
    if(chain %in% c("TCRA", "TCRG")) {
        ref <- 1
    } else if(chain %in% c("TCRB", "TCRD")) {
        ref <- 2
    }
    ref2 <- paste0("cdr3_", sequence, ref)
    bound <- bind_rows(df)
    #Should make it work as either grouped or non-grouped
    if (!is.null(group)) {
      bound <- split(bound, bound[,group])
      list.length <- length(bound)
    } else {
      bound <- list(bound)
      list.length <- 1
    }
    
    for (x in seq_along(bound)) {
      dictionary <- na.omit(unique(bound[[x]][,ref2]))
      dictionary <- str_split(dictionary, ";", simplify = TRUE)[,1]
      length <- nchar(dictionary)
      matrix <- as.matrix(stringdistmatrix(dictionary, method = "lv"))    
      out_matrix <- matrix(ncol = ncol(matrix), nrow=ncol(matrix))
      for (j in seq_len(ncol(matrix))) {
          for (k in seq_len(nrow(matrix))) {
              if (j == k) {
                  out_matrix[j,k] <- NA
              } else{
                  if (length[j] - length[k] >= round(mean(length)/2)) {
                      out_matrix[j,k] <- matrix[j,k]/(max(length[j], length[k]))
                      out_matrix[k,j] <- matrix[k,j]/(max(length[j], length[k]))
                  } else {
                  out_matrix[j,k] <- matrix[j,k]/((length[j]+ length[k])/2)
                  out_matrix[k,j] <- matrix[k,j]/((length[j]+ length[k])/2)
                  }
              }
          }
      }
      filtered <- which(out_matrix <= (1-threshold), arr.ind = TRUE)
      if (nrow(filtered) > 0) { 
          for (i in seq_len(nrow(filtered))) {
              max <- max(filtered[i,])
              min <- min(filtered[i,])
              filtered[i,1] <- max
              filtered[i,2] <- min
          }
          filtered <- unique(filtered) #removing redundant comparisons
          out <- NULL
          colnames(filtered) <- c("To", "From")
          
          g <- graph_from_data_frame(filtered)
          components <- components(g, mode = c("weak"))
          out <- data.frame("cluster" = components$membership, 
                            "filtered" = names(components$membership))
          if (list.length == 1) {
            out$cluster <- paste0(chain, ":LD", ".", out$cluster)
          } else {
            out$cluster <- paste0(names(bound)[x], ".", chain, ":LD", ".", out$cluster)
          }
          out$filtered <- dictionary[as.numeric(out$filtered)]
          
          uni_IG <- as.data.frame(unique(dictionary[dictionary %!in% out$filtered]))
          colnames(uni_IG) <- "filtered"
          if (list.length == 1) {
            uni_IG$cluster <- paste0(chain, ".", seq_len(nrow(uni_IG))) 
          } else {
            uni_IG$cluster <- paste0(names(bound)[x], ".", chain, ".", seq_len(nrow(uni_IG))) 
          }
      }
      output <- rbind.data.frame(out, uni_IG)
      colname <- paste0(chain, "_cluster")
      colnames(output) <- c(colname, ref2)
      output.list[[x]] <- output
    }
      for (i in seq_along(df)) {
          tmp <- df[[i]]
          output <- bind_rows(output.list)
          tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = TRUE)[,1]
          output2 <- output[output[,2] %in% tmp[,ref2],]
          
          tmp <-  suppressMessages(join(tmp,  output2))
          df[[i]] <- tmp
      }
    return(df)
}  

