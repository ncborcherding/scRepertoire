#' Examining the diversity of amino acids by position
#'
#' This function the diversity amino acids along the residues 
#' of the CDR3 amino acid sequence. Please see 
#' [clonalDiversity()] for more information on 
#' the underlying methods for diversity/entropy calculations. 
#' Positions without variance will have a value reported as 0 
#' for the purposes of comparison.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' positionalEntropy(combined, 
#'                   chain = "TRB", 
#'                   aa.length = 20)

#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param chain indicate a specific chain should be used - 
#' e.g."TRA", "TRG", "IGH", "IGL", etc
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order for `group.by` or 
#' "alphanumeric" to plot groups in order
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method The method to calculate the entropy/diversity - 
#' `"shannon"`, `"inv.simpson"`, `"gini.simpson"`, `"norm.entropy"`, 
#' `"pielou"`, `"hill0"`, `"hill1"`, `"hill2"`
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any [hcl.pals][grDevices::hcl.pals]
#'
#' @export
#' @importFrom immApex calculateEntropy
#' @concept Summarize_Repertoire
#' @return ggplot of line graph of diversity by amino acid residue position
positionalEntropy <- function(input.data, 
                              chain = "TRB", 
                              group.by = NULL, 
                              order.by = NULL,
                              aa.length = 20,
                              method = "norm.entropy",
                              exportTable = FALSE, 
                              palette = "inferno")  {
  
  sco <- .is.seurat.or.se.object(input.data)
  input.data <- .dataWrangle(input.data, 
                              group.by, 
                              .theCall(input.data, "CTaa", check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, "CTaa")
  
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  #Selecting Diversity Function
  lapply(seq_along(input.data), function(x){
    strings <- .processStrings(input.data[[x]][["CTaa"]], aa.length)
    entropy <- calculateEntropy(strings, method = method)
    entropy_df <- as.data.frame(entropy)
    entropy_df$Position <- 1:nrow(entropy_df)
    entropy_df$group <- names(input.data)[x]
    rownames(entropy_df) <- NULL # Remove the messy row names
    entropy_df
  }) -> group.results
  
  mat_melt <- do.call(rbind, group.results)
  
  if(!is.null(order.by)) {
    mat_melt <- .orderingFunction(vector = order.by,
                                  group.by = "group", 
                                  mat_melt)
  }
    
  plot <- ggplot(mat_melt, aes(x=as.factor(Position), 
                               y = .data[["entropy"]], 
                               group= group, 
                               color = group)) +
          geom_line(stat = "identity") +
          geom_point() + 
          scale_color_manual(name = "Groups", 
                            values = rev(.colorizer(palette,length(input.data)))) +
          xlab("Amino Acid Residues") +
          ylab("Relative Diversity") +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if (exportTable == TRUE) { 
      return(mat_melt) 
    }
    return(plot)
}

  


