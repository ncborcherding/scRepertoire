#' Examining the relative amino acid composition by position
#'
#' This function the proportion of amino acids along the residues 
#' of the CDR3 amino acid sequence.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentAA(combined, 
#'           chain = "TRB", 
#'           aa.length = 20)
  
#' @param input.data The product of [combineTCR()], [combineBCR()], or
#'  [combineExpression()].
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL".
#' @param group.by The variable to use for grouping.
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any [hcl.pals][grDevices::hcl.pals].
#' @importFrom dplyr mutate_at mutate_if
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of stacked bar graphs of amino acid proportions
percentAA <- function(input.data, 
                      chain = "TRB", 
                      group.by = NULL, 
                      order.by = NULL,
                      aa.length = 20,
                      exportTable = FALSE, 
                      palette = "inferno")  {
  
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
  input.data <- .data.wrangle(input.data, group.by, "CTaa", chain)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  #Getting AA Counts
  aa.count.list <- .aa.counter(input.data, "CTaa", aa.length)
  
  #Calculating proportion and melting data
  lapply(seq_along(aa.count.list), function(x) {
    aa.count.list[[x]] <- aa.count.list[[x]] %>% mutate_if(is.numeric, list(~ ./sum(.)))
    melt.res <- suppressMessages(melt(aa.count.list[[x]]))
    melt.res$group <- names(input.data)[x]
    melt.res
  }) -> res.list
  
  mat_melt <- do.call(rbind, res.list)
  if(!is.null(order.by)) {
    mat_melt <- .ordering.function(vector = order.by,
                                   group.by = "variable", 
                                   mat_melt)
  }
  
  plot <- ggplot(mat_melt, aes(x=as.factor(variable), y = value, fill=AA)) +
    geom_bar(stat = "identity", position="fill", lwd= 0.25, color = "black") +
    scale_fill_manual(name = "Amino Acid", 
                      values = rev(.colorizer(palette,21))) +
    xlab("Amino Acid Residues") +
    ylab("Relative Percent") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if(length(res.list) > 1) {
    plot <- plot + facet_grid(group~.)
  }
  if (exportTable == TRUE) { 
    return(mat_melt) 
  }
  return(plot)
}    
    
