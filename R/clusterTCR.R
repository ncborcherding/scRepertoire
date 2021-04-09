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
#' sub_combined <- clusterTCR(combined[[1]], chain = "TCRA", sequence = "aa")
#' 
#' @param df The product of CombineTCR() or CombineBCR().
#' @param chain The TCR to cluster
#' @param sequence Clustering based on either "aa" or "nt"
#' @param threshold The normalized edit distance to consider. The higher the number the more 
#' similarity of sequence will be used for clustering.
#' @importFrom stringdist stringdistmatrix
#' @importFrom igraph graph_from_data_frame components
#' @importFrom plyr join
#' @export
#' @return List of clonotypes for individual cell barcodes

clusterTCR <- function(df, chain = NULL, sequence = NULL, threshold = 0.85) {
  `%!in%` = Negate(`%in%`)
  df <- checkList(df)
  if(chain %in% c("TCRA", "TCRG")) {
    ref <- 1
  } else if(chain %in% c("TCRB", "TCRD")) {
    ref <- 2
  }
  ref2 <- paste0("cdr3_", sequence, ref)
  bound <- dplyr::bind_rows(df)
  dictionary <- na.omit(unique(bound[,ref2]))
  dictionary <- str_split(dictionary, ";", simplify = T)[,1]
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
        }
        out_matrix[j,k] <- matrix[j,k]/((length[j]+ length[k])/2)
        out_matrix[k,j] <- matrix[k,j]/((length[j]+ length[k])/2)
      }
    }
  }
  filtered <- which(out_matrix <= (1-threshold), arr.ind = TRUE)
  if (nrow(filtered) > 0) { 
    for (i in 1:nrow(filtered)) {
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
    out$cluster <- paste0(chain, ":LD", ".", out$cluster)
    out$filtered <- dictionary[as.numeric(out$filtered)]
    
    uni_IG <- as.data.frame(unique(dictionary[dictionary %!in% out$filtered]))
    colnames(uni_IG) <- "filtered"
    uni_IG$cluster <- paste0(chain, ".", seq_len(nrow(uni_IG)))
  }
  output <- rbind.data.frame(out, uni_IG)
  colname <- paste0(chain, "_cluster")
  colnames(output) <- c(colname, ref2)
  for (i in seq_along(df)) {
    tmp <- df[[i]]
    tmp[,ref2] <- str_split(tmp[,ref2], ";", simplify = T)[,1]
    output2 <- output[output[,2] %in% tmp[,ref2],]
    
    tmp <-  suppressMessages(join(tmp,  output2))
    df[[i]] <- tmp
  }
  return(df)
} 