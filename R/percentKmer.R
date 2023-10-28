#' Examining the relative composition of kmer motifs
#'
#' This function the of kmer for nucleotide or amino acids 
#' in the CDR3 sequence.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentKmer(combined, chain = "TRB", motif.length = 3)

#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL".
#' @param cloneCall How to call the clonotype - CDR3 nucleotide (nt) or 
#' CDR3 amino acid (aa).
#' @param group.by The variable to use for grouping.
#' @param motif.length The length of the kmer to analyze.
#' @param top.motifs Return the n most variable motifs as a function of 
#' median absolute deviation. 
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats mad
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of percentage of kmers as a heatmap

percentKmer <- function(df, 
                        chain = "TRB", 
                        cloneCall = "aa",
                        group.by = NULL, 
                        motif.length = 3,
                        top.motifs = 30,
                        exportTable = FALSE, 
                        palette = "inferno") {
  motifs.to.save <- NULL
  cloneCall <- .theCall(cloneCall)
  sco <- is_seurat_object(df) | is_se_object(df)
  df <- .data.wrangle(df, group.by, cloneCall, chain)
  if(!is.null(group.by) && !sco) {
    df <- .groupList(df, group.by)
  }
  #Determining the function to generate motifs
  kmerFunc <- switch(cloneCall,
                           "CTaa" = .generate_unique_aa_motifs,
                           "CTnt" = .generate_unique_nt_motifs,
                           stop("Please select 'nt' or 'aa' for the cloneCall"))
  #Getting motifs and matrix to place the counts
  unique.motifs <- kmerFunc(motif.length)
  mat <- matrix(NA, ncol = length(unique.motifs), nrow = length(df))
  colnames(mat) <- unique.motifs
  rownames(mat) <- names(df)
  
  #Looping through the counts of motifs
  for (i in seq_along(df)) {
    df[[i]][,cloneCall][which(df[[i]][, cloneCall] == "NA")] <- NA
    df[[i]] <- df[[i]][!is.na(df[[i]][, cloneCall]),]

    # if (identical(cloneCall, "CTnt")) {
    #   mat[i, ] <- rcppGetNtKmerPercent(df[[i]][, cloneCall], motif.length)
    #   next
    # }

    motifs <- .tokenize_multiple_sequences(df[[i]][,cloneCall], motif.length)
    motif.table <- table(unlist(motifs))
    if(any(grepl("\\;", names(motif.table)))) {
      motif.table <- motif.table[!grepl("\\;", names(motif.table))]
    }
    mat.pos <- match(names(motif.table), colnames(mat))
    mat[i, mat.pos] <- as.vector(motif.table)
    mat[i, ] <- mat[i, ] / sum(na.omit(mat[i, ]))
  }

  #Removing any column that is all NAs
  if(any(colSums(is.na(mat)) == length(df))) {
    mat <- mat[,-which(colSums(is.na(mat)) == length(df))]
  }
  mat[is.na(mat)] <- 0
  
  #Filtering for top.motifs
  if(!is.null(top.motifs)) {
    mads <- apply(mat, 2, mad)
    motifs.to.save <- names(sort(mads, decreasing = TRUE))[seq_len(top.motifs)]
    mat <- mat[, colnames(mat) %in% motifs.to.save]
  }

  #Getting mat into a ggplot-compliant form
  mat_melt <- melt(mat)
  if (!is.null(motifs.to.save)) {
    mat_melt$Var2 <- factor(mat_melt$Var2, levels = rev(motifs.to.save))
  }
  
  if (exportTable) {
    return(mat)
  }

  #Plotting
  ggplot(mat_melt, aes(x=Var2, y = Var1, fill=value)) +
    geom_tile(lwd= 0.1, color = "black") + 
    scale_fill_gradientn(name = "Percentage", colors = .colorizer(palette,21)) +
    theme_classic() + 
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_blank())
}

.generate_unique_aa_motifs <- function(motif_length) {
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  all_motifs <- expand.grid(replicate(motif_length, amino_acids, simplify = FALSE))
  unique_motifs <- unique(apply(all_motifs, 1, paste, collapse = ""))
  return(unique_motifs)
}

.generate_unique_nt_motifs <- function(motif_length) {
  rcppGenerateUniqueNtMotifs(motif_length)
}

.tokenize_sequence <- function(sequence, motif_length) {
  tokens <- c()
  sequence_length <- nchar(sequence)
  
  for (i in 1:(sequence_length - motif_length + 1)) {
    motif <- substr(sequence, i, i + motif_length - 1)
    tokens <- c(tokens, motif)
  }
  
  return(tokens)
}

.tokenize_multiple_sequences <- function(sequences, motif_length) {
  sapply(sequences, .tokenize_sequence, motif_length = motif_length)
}
