#' Examining the V or J gene usage across groupings
#'
#' This function the proportion V or J genes used by 
#' grouping variables
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentGenes(combined, chain = "TRB", gene = "Vgene)

#' @param df The product of combineTCR(), combineBCR(), 
#' expression2List(), or combineExpression().
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param gene "V", "D" or "J"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2
#' @importFrom stringr str_split str_sort 
#' @importFrom reshape2 melt
#' @export
#' @return ggplot of percentage of indicated genes as a heatmap
#' 
percentGenes <- function(df,
                         chain = "TRB",
                         group.by = NULL, 
                         split.by = NULL,
                         exportTable = FALSE, 
                         palette = "inferno") {
  
  df <- list.input.return(df, split.by)
  df <- checkBlanks(df, "CTgene")
  for(i in seq_along(df)) {
    df[[i]] <- off.the.chain(df[[i]], chain, "CTgene")
  }
  
  if(chain %in% c("TRA", "TRG", "IGL")) {
    positions <- c(1,2)
  } else {
    positions <- c(1,3)
  }
  
  
  #Getting indicated genes across list
  gene_counts <- lapply(df, function(x) {
      tmp <- unlist(str_split(x[,"CTgene"], ";"))
      tmp <- str_split(tmp, "[.]", simplify = TRUE)
      strings <- paste0(tmp[,positions[1]], ";", tmp[,positions[2]])
      strings <- strings[-which(strings == "NA;")]
  })
  #Need total unique genes
  gene.dictionary <- unique(unlist(gene_counts))
  
  #Summarizing the gene usage across the list
  summary <- lapply(gene_counts, function(x) {
                 percentages <- unlist(prop.table(table(x)))
                 genes.to.add <- gene.dictionary [which(gene.dictionary  %!in% names(percentages))]
                 if(length(genes.to.add) >= 1) {
                   percentages.to.add <- rep(0, length(genes.to.add))
                   names(percentages.to.add) <- genes.to.add
                   percentages <- c(percentages, percentages.to.add)
                 }
                 percentages <- percentages[match(str_sort(names(percentages), numeric = TRUE), names(percentages))]
                 
  
  })
  summary <- do.call(rbind,summary)
  if (exportTable == TRUE) { 
    return(summary) 
  }
  #Melting matrix and Visualizing
  mat_melt <- melt(summary)
  plot <- ggplot(mat_melt, aes(y=Var1, x = Var2, fill=round(value*100,2))) +
    geom_tile(lwd= 0.25, color = "black") +
    scale_fill_gradientn(name = "Percentage", colors = .colorizer(palette,21)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.title = element_blank())
  return(plot)
}
