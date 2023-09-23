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
#' combined <- combineTCR(contig_list, 
#'                        samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                    "P19B","P19L", "P20B", "P20L"))
#' clonalCompare(combined, 
#'               top.clones = 5, 
#'               samples = c("P17B", "P17L"), 
#'               cloneCall="aa")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param samples The specific samples to isolate for visualization.
#' @param clonotypes The specific clonal sequences of interest.
#' @param top.clones The top number of clonotype sequences per group
#' @param highlight.clones Clonal sequences to highlight, if present, 
#' all other clones returned will be grey.
#' @param relabel.clones Simplify the legend of the graph by returning
#' clones that are numerically indexed.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param graph The type of graph produced, either "alluvial" or "area".
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2
#'
#' @export
#' @return ggplot of the proportion of total sequencing read of 
#' selecting clonotypes
clonalCompare <- function(df, 
                              cloneCall = "strict", 
                              chain = "both", 
                              samples = NULL, 
                              clonotypes = NULL, 
                              top.clones = NULL,
                              highlight.clones = NULL,
                              relabel.clones = FALSE,
                              split.by = NULL, 
                              graph = "alluvial", 
                              exportTable = FALSE, 
                              palette = "inferno"){
  
  #Tie goes to indicated clones over top clones
  if(!is.null(top.clones) & !is.null(clonotypes)) {
    top.clones <- NULL
  }
  
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  Con.df <- NULL
  
  #Loop through the list to get a proportional summary
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
  
  #Filtering steps 
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,] 
  }
  if (!is.null(clonotypes)) {
    Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] 
  } else if (!is.null(top.clones)) {
    top <- Con.df %>%
      group_by(Con.df[,3]) %>%
      slice_max(n = top.clones, order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] 
  }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") 
  }
  #Clonotype relabeling
  clones.returned <- as.vector(unique(Con.df[,"Clonotypes"]))
  if (relabel.clones) {
    new.clones <- paste0("Clonotype: ", seq_len(length(clones.returned)))
    names(new.clones) <- clones.returned
    #Isolated new clone names for highlight purposes
    if(!is.null(highlight.clones)) {
      highlight.clones <- unname(new.clones[which(names(new.clones) %in% highlight.clones)])
    }
    Con.df[,"Clonotypes"] <- new.clones[as.vector(Con.df[,"Clonotypes"])]
    clones.returned <- as.vector(unique(Con.df[,"Clonotypes"]))
  }
  if (exportTable == TRUE) { 
    return(Con.df)
  }
  
  #Plotting Functions
  plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                             stratum = Clonotypes, alluvium = Clonotypes, 
                             y = Proportion, label = Clonotypes)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), 
          legend.text=element_text(size=rel(0.5)), 
          legend.key.size = unit(0.5,"line"))
  if (graph == "alluvial") {
    plot <- plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot <- plot +
      geom_area(aes(group = Clonotypes), color = "black")
  }
  
  #Highlighting specific clones
  if (!is.null(highlight.clones)) {
    clone.colors <- rep("grey", length(clones.returned))
    pos <- which(clones.returned %in% highlight.clones)
    clone.colors[pos] <- .colorizer(palette, length(pos))
    names(clone.colors) <- clones.returned
    plot <- plot + scale_fill_manual(values = clone.colors)
  } else {
    plot <- plot + scale_fill_manual(values = .colorizer(palette, length(unique(Con.df[,"Clonotypes"]))))
  }
  return(plot)
}
