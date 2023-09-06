#' Hierarchical clustering of clonotypes on clonotype size and 
#' Jensen-Shannon divergence
#'
#' This function produces a hierarchical clustering of clonotypes by sample 
#' using the Jensen-Shannon distance and discrete gamma-GPD spliced threshold 
#' model in 
#' \href{https://bioconductor.org/packages/devel/bioc/html/powerTCR.html}{powerTCR R package}.
#' Please read and cite PMID: 30485278 if using the function for analyses. 
#' 
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' clonesizeDistribution(combined, cloneCall = "strict", method="ward.D2")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param threshold Numerical vector containing the thresholds 
#' the grid search was performed over.
#' @param method The clustering parameter for the dendrogram.
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph.
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot
#' @importFrom powerTCR fdiscgammagpd get_distances
#' @export
#' @return ggplot dendrogram of the clone size distribution

clonesizeDistribution <- function(
    df,
    cloneCall ="strict", 
    chain = "both", 
    method = "ward.D2", 
    threshold = 1, 
    group.by = NULL,
    split.by = NULL, 
    exportTable = FALSE
) {
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  if(!is.null(group.by)) {
    df <- groupList(df, group.by)
  }
  data <- bind_rows(df)
  unique_df <- unique(data[,cloneCall])
  Con.df <- data.frame(matrix(NA, length(unique_df), length(df)))
  Con.df <- data.frame(unique_df, Con.df, stringsAsFactors = FALSE)
  colnames(Con.df)[1] <- "clonotype"
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
    data <- df[[i]]
    data <- data.frame(table(data[,cloneCall]), 
                       stringsAsFactors = FALSE)
    colnames(data) <- c(cloneCall, "Freq")
    for (y in seq_along(unique_df)){
      clonotype.y <- Con.df$clonotype[y]
      location.y <- which(clonotype.y == data[,cloneCall])
      Con.df[y,i+1] <- data[location.y[1],"Freq"]
    }
  }
  colnames(Con.df)[2:(length(df)+1)] <- names(df)
  Con.df[is.na(Con.df)] <- 0
  list <- list()
  for (i in seq_along(df)) {
    list[[i]] <- Con.df[,i+1]
    list[[i]] <- suppressWarnings(fdiscgammagpd(list[[i]], useq = threshold))
  }
  names(list) <- names(df)
  grid <- 0:10000
  distances <- get_distances(list, grid, modelType="Spliced")
  hclust <- hclust(as.dist(distances), method = method)
  hcd <- as.dendrogram(hclust)
  plot <- plot(hcd)
  if (exportTable) { return(distances) }
  return(plot)
}