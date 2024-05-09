#' Examining the diversity of amino acids by position
#'
#' This function the diversity amino acids along the residues 
#' of the CDR3 amino acid sequence. Please see 
#' \code{\link{clonalDiversity}} for more information on 
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

#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method The method to calculate the entropy/diversity - 
#' "shannon", "inv.simpson", "norm.entropy"
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}
#' @import ggplot2
#' @importFrom stringr str_split
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of line graph of diversity by position
positionalEntropy <- function(input.data, 
                              chain = "TRB", 
                              group.by = NULL, 
                              order.by = NULL,
                              aa.length = 20,
                              method = "norm.entropy",
                              exportTable = FALSE, 
                              palette = "inferno")  {
  
  if(method %!in% c("shannon", "inv.simpson", "norm.entropy")) {
    stop("Please select a compatible method.")
  }
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, "CTaa", check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, "CTaa")
  
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  #Selecting Diversity Function
  diversityFunc <- switch(method,
                          "norm.entropy" = .normentropy,
                          "inv.simpson" = .invsimpson,
                          "shannon" = .shannon,
                          stop("Invalid method provided"))
  
  aa.count.list <- .aa.counter(input.data, "CTaa", aa.length)
  
  lapply(aa.count.list, function(x){
      diversity <- sapply(x[,2:ncol(x)], diversityFunc)
      diversity[is.nan(diversity)] <- 0
      diversity
  }) -> group.results

  mat <- do.call(rbind, group.results)
  mat_melt <- suppressMessages(melt(mat))
  
  if(!is.null(order.by)) {
    mat_melt <- .ordering.function(vector = order.by,
                                   group.by = "Var1", 
                                   mat_melt)
  }
    
  plot <- ggplot(mat_melt, aes(x=Var2, y = value, group= Var1, color = Var1)) +
          geom_line(stat = "identity") +
          geom_point() + 
          scale_color_manual(name = "Groups", 
                            values = rev(.colorizer(palette,nrow(mat)))) +
          xlab("Amino Acid Residues") +
          ylab("Relative Diversity") +
          theme_classic() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    if (exportTable == TRUE) { 
      return(mat_melt) 
    }
    return(plot)
}

  


