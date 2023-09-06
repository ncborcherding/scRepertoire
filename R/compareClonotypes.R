#' Demonstrate the difference in clonal proportion between clonotypes
#'
#' This function produces an alluvial or area graph of the proportion of 
#' the indicated clonotypes for all or selected samples. Clonotypes can be 
#' selected using the clonotypes parameter with the specific sequence of 
#' interest or using the number parameter with the top n clonotypes by 
#' proportion to be visualized. 
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' compareClonotypes(combined, numbers = 10, 
#' samples = c("PX_P", "PX_T"), cloneCall="aa")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param samples The specific samples to isolate for visualization.
#' @param clonotypes The specific sequences of interest.
#' @param numbers The top number clonotype sequences per group
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param graph The type of graph produced, either "alluvial" or "area".
#' @param exportTable Returns the data frame used for forming the graph.
#' @import ggplot2
#'
#' @export
#' @return ggplot of the proportion of total sequencing read of 
#' selecting clonotypes
compareClonotypes <- function(df, 
                              cloneCall = "strict", 
                              chain = "both", 
                              samples = NULL, 
                              clonotypes = NULL, 
                              numbers = NULL, 
                              split.by = NULL, 
                              graph = "alluvial", 
                              exportTable = FALSE){
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  if (!is.null(numbers) & !is.null(clonotypes)) {
    stop("Make sure your inputs are either numbers or clonotype sequences.")
  }
  Con.df <- NULL
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
    tbl <- as.data.frame(table(df[[i]][,cloneCall]))
    tbl[,2] <- tbl[,2]/sum(tbl[,2])
    colnames(tbl) <- c("Clonotypes", "Proportion")
    tbl$Sample <- names(df[i])
    Con.df <- rbind.data.frame(Con.df, tbl)
  }
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,] }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
  if (!is.null(numbers)) {
    top <- Con.df %>%
      group_by(Con.df[,3]) %>%
      slice_max(n = numbers, order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") }
  if (exportTable == TRUE) { return(Con.df)}
  
  plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                             stratum = Clonotypes, alluvium = Clonotypes, 
                             y = Proportion, label = Clonotypes)) +
    theme_classic() +
    theme(axis.title.x = element_blank())
  if (graph == "alluvial") {
    plot <- plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot <- plot +
      geom_area(aes(group = Clonotypes), color = "black") }
  return(plot)
}