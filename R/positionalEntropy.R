#' Examining the Diversity of Amino Acids by Position
#'
#' This function the diversity amino acids along the residues of the CDR3 
#' amino acid sequence. Please see [clonalDiversity()] for more information 
#' on the underlying methods for diversity/entropy calculations. Positions 
#' without variance will have a value reported as 0 for the purposes of comparison.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using positionalEntropy()
#' positionalEntropy(combined, 
#'                   chain = "TRB", 
#'                   aa.length = 20)
#'                   
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param chain The TCR/BCR chain to use. Accepted values: `TRA`, `TRB`, `TRG`, 
#' `TRD`, `IGH`, or `IGL` (for both light chains).
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method The method to calculate the entropy/diversity - 
#' `"shannon"`, `"inv.simpson"`, `"gini.simpson"`, `"norm.entropy"`, 
#' `"pielou"`, `"hill0"`, `"hill1"`, `"hill2"`
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#' 
#' @export
#' @importFrom immApex calculateEntropy
#' @concept Summarize_Repertoire
#' @return A ggplot object displaying entropy or diversity by amino acid position.
#' If `exportTable = TRUE`, a matrix of the raw data is returned.
positionalEntropy <- function(input.data, 
                              chain = "TRB", 
                              group.by = NULL, 
                              order.by = NULL,
                              aa.length = 20,
                              method = "norm.entropy",
                              exportTable = FALSE, 
                              palette = "inferno",
                              ...)  {
  
  sco <- .is.seurat.or.se.object(input.data)
  input.data <- .dataWrangle(input.data, 
                              group.by, 
                              .theCall(input.data, "CTaa", 
                                       check.df = FALSE, silent = TRUE), 
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
                               group= .data[["group"]], 
                               color = group)) +
          geom_line(stat = "identity") +
          geom_point() + 
          scale_color_manual(name = "Groups", 
                            values = rev(.colorizer(palette,length(input.data)))) +
          xlab("Amino Acid Residues") +
          ylab("Relative Diversity") +
          .themeRepertoire(...) + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if (exportTable == TRUE) { 
      return(mat_melt) 
    }
    return(plot)
}

  


