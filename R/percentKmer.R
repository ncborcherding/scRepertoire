#' Examining the relative composition of kmer motifs in clones.
#'
#' This function the of kmer for nucleotide (**nt**) or 
#' amino acid (**aa**) sequences. Select the length of the
#' kmer to quantify using the **motif.length** parameter.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentKmer(combined, 
#'             chain = "TRB", 
#'             motif.length = 3)
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param cloneCall How to call the clone - CDR3 nucleotide (**nt**) or 
#' CDR3 amino acid (**aa**)
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param motif.length The length of the kmer to analyze
#' @param top.motifs Return the n most variable motifs as a function of 
#' median absolute deviation
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @importFrom stats mad
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of percentage of kmers as a heatmap

percentKmer <- function(input.data, 
                        chain = "TRB", 
                        cloneCall = "aa",
                        group.by = NULL, 
                        order.by = NULL,
                        motif.length = 3,
                        top.motifs = 30,
                        exportTable = FALSE, 
                        palette = "inferno") {
  
  if(cloneCall %!in% c("aa", "nt")) {
    stop("Please select either nucleotide (nt) or amino acid (aa) sequences for cloneCall")
  }
  motifs.to.save <- NULL
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
  if(!is.null(group.by) && !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  #Determining the function to generate motifs
  kmerFunc <- switch(cloneCall,
                           "CTaa" = .generate_unique_aa_motifs,
                           "CTnt" = .generate_unique_nt_motifs,
                           stop("Please select 'nt' or 'aa' for the cloneCall"))
  #Getting motifs and matrix to place the counts
  unique.motifs <- kmerFunc(motif.length)
  mat <- matrix(NA, ncol = length(unique.motifs), nrow = length(input.data))
  colnames(mat) <- unique.motifs
  rownames(mat) <- names(input.data)
  
  #Looping through the counts of motifs
  for(i in seq_along(input.data)) {
    input.data[[i]][,cloneCall][which(input.data[[i]][,cloneCall] == "NA")] <- NA
    input.data[[i]] <- input.data[[i]][!is.na(input.data[[i]][,cloneCall]),]

    if (identical(cloneCall, "CTnt") && motif.length > 1L && motif.length <= 32L) {
      mat[i, ] <- rcppGetNtKmerPercent(input.data[[i]][,cloneCall], motif.length)
      next
    }

    if (identical(cloneCall, "CTaa") && motif.length > 1L && motif.length <= 12L) {
      mat[i, ] <- rcppGetAaKmerPercent(input.data[[i]][,cloneCall], unique.motifs, motif.length)
      next
    }

    # this point could only be reached if motif.length is extremely long or is one.
    # which wouldn't really be practical. Maybe the following should just be deleted
    # and a cap could be put on motif.length?

    motifs <- .tokenize_multiple_sequences(input.data[[i]][,cloneCall], motif.length)
    motif.table <- table(unlist(motifs))
    if(any(grepl("\\;", names(motif.table)))) {
      motif.table <- motif.table[!grepl("\\;", names(motif.table))]
    }
    mat.pos <- match(names(motif.table), colnames(mat))
    mat[i, mat.pos] <- as.vector(motif.table)
    mat[i, ] <- mat[i, ] / sum(na.omit(mat[i, ]))
  }

  #Removing any column that is all NAs
  if(any(colSums(is.na(mat)) == length(input.data))) {
    mat <- mat[,-which(colSums(is.na(mat)) == length(input.data))]
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
  
  if(!is.null(order.by)) {
    mat_melt <- .ordering.function(vector = order.by,
                                   group.by = "Var1", 
                                   mat_melt)
  }
  
  if (exportTable) {
    return(mat)
  }

  #Plotting
  plot <- ggplot(mat_melt, aes(x=Var2, y = Var1, fill=value)) +
            geom_tile(lwd = 0.1, color = "black") + 
            scale_fill_gradientn(name = "Percentage", colors = .colorizer(palette,21)) +
            theme_classic() + 
            coord_flip() + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
                  axis.title = element_blank())
  return(plot)
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
  
  for (i in seq_len((sequence_length - motif_length + 1))) {
    motif <- substr(sequence, i, i + motif_length - 1)
    tokens <- c(tokens, motif)
  }
  
  return(tokens)
}

.tokenize_multiple_sequences <- function(sequences, motif_length) {
  sapply(sequences, .tokenize_sequence, motif_length = motif_length)
}
