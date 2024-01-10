#' Demonstrate the difference in clonal proportion between clones
#'
#' This function produces an alluvial or area graph of the proportion of 
#' the indicated clones for all or selected samples (using the 
#' \strong{samples} parameter). Individual clones can be selected 
#' using the \strong{clones} parameter with the specific sequence of 
#' interest or using the \strong{top.clones} parameter with the top 
#' n clones by proportion to be visualized. 
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
#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}.
#' @param cloneCall How to call the clone - VDJC gene (\strong{gene}), 
#' CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
#' VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable 
#' in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param samples The specific samples to isolate for visualization.
#' @param clones The specific clonal sequences of interest.
#' @param top.clones The top number of clonal sequences per group.
#' (e.g., top.clones = 5)
#' @param highlight.clones Clonal sequences to highlight, if present, 
#' all other clones returned will be grey.
#' @param relabel.clones Simplify the legend of the graph by returning
#' clones that are numerically indexed.
#' @param group.by If using a single-cell object, the column header
#' to group the new list. \strong{NULL} will return the active 
#' identity or cluster.
#' @param graph The type of graph produced, either \strong{"alluvial"} 
#' or \strong{"area"}.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom stringr str_sort
#'
#' @export
#' @concept Visualizing_Clones
#' @return ggplot of the proportion of total sequencing read of 
#' selecting clones
clonalCompare <- function(input.data, 
                          cloneCall = "strict", 
                          chain = "both", 
                          samples = NULL, 
                          clones = NULL, 
                          top.clones = NULL,
                          highlight.clones = NULL,
                          relabel.clones = FALSE,
                          group.by = NULL, 
                          graph = "alluvial", 
                          exportTable = FALSE, 
                          palette = "inferno"){
  
  #Tie goes to indicated clones over top clones
  if(!is.null(top.clones) & !is.null(clones)) {
    top.clones <- NULL
  }
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
  
  Con.df <- NULL
  
  #Loop through the list to get a proportional summary
  for (i in seq_along(input.data)) {
    if (chain != "both") {
      input.data[[i]] <- .off.the.chain(input.data[[i]], chain, cloneCall)
    }
    tbl <- as.data.frame(table(input.data[[i]][,cloneCall]))
    tbl[,2] <- tbl[,2]/sum(tbl[,2])
    colnames(tbl) <- c("clones", "Proportion")
    tbl$Sample <- names(input.data[i])
    Con.df <- rbind.data.frame(Con.df, tbl)
  }
  
  #Filtering steps 
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,] 
  }
  if (!is.null(clones)) {
    Con.df <- Con.df[Con.df$clones %in% clones,] 
  } else if (!is.null(top.clones)) {
    top <- Con.df %>%
      group_by(Con.df[,3]) %>%
      slice_max(n = top.clones, order_by = Proportion, with_ties = FALSE)
    Con.df <- Con.df[Con.df$clones %in% top$clones,] 
  }
  if (nrow(Con.df) < length(unique(Con.df$Sample))) {
    stop("Reasses the filtering strategies here, there are not 
            enough clones to examine.") 
  }
  #Clones relabeling
  clones.returned <- as.vector(unique(Con.df[order(Con.df[,"Proportion"], decreasing = TRUE),"clones"]))
  if (relabel.clones) {
    new.clones <- paste0("Clone: ", seq_len(length(clones.returned)))
    names(new.clones) <- clones.returned
    #Isolated new clone names for highlight purposes
    if(!is.null(highlight.clones)) {
      highlight.clones <- unname(new.clones[which(names(new.clones) %in% highlight.clones)])
    }
    Con.df[,"clones"] <- new.clones[as.vector(Con.df[,"clones"])]
    Con.df[,"clones"] <- factor(Con.df[,"clones"], 
                                levels = str_sort(unique(Con.df[,"clones"]), numeric = TRUE))
    clones.returned <- as.vector(unique(Con.df[,"clones"]))
  }
  if (exportTable == TRUE) { 
    return(Con.df)
  }
  
  #Plotting Functions
  plot <- ggplot(Con.df, aes(x = Sample, 
                             fill = clones, 
                             group = clones,
                             stratum = clones, 
                             alluvium = clones, 
                             y = Proportion, 
                             label = clones)) +
    theme_classic() +
    theme(axis.title.x = element_blank(), 
          legend.text=element_text(size=rel(0.5)), 
          legend.key.size = unit(0.5,"line"))
  if (graph == "alluvial") {
    plot <- plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot <- plot +
      geom_area(aes(group = clones), color = "black")
  }
  
  #Highlighting specific clones
  if (!is.null(highlight.clones)) {
    clone.colors <- rep("grey", length(clones.returned))
    pos <- which(clones.returned %in% highlight.clones)
    clone.colors[pos] <- .colorizer(palette, length(pos))
    names(clone.colors) <- clones.returned
    plot <- plot + scale_fill_manual(values = clone.colors)
  } else {
    plot <- plot + scale_fill_manual(values = .colorizer(palette, length(unique(Con.df[,"clones"]))))
  }
  return(plot)
}
