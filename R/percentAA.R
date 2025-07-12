#' Examining the relative amino acid composition by position
#'
#' This function the proportion of amino acids along the residues 
#' of the CDR3 amino acid sequence.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentAA(combined, 
#'           chain = "TRB", 
#'           aa.length = 20)
  
#' @param input.data The product of [combineTCR()], [combineBCR()], or
#'  [combineExpression()].
#' @param chain indicate a specific chain should be used - 
#' e.g. "TRA", "TRG", "IGH", "IGL", etc
#' @param group.by The variable to use for grouping.
#' @param order.by A vector of specific plotting order for `group.by` or 
#' "alphanumeric" to plot groups in order
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any [hcl.pals][grDevices::hcl.pals].
#' @param ... Additional arguments passed to the ggplot theme
#' @importFrom immApex calculateFrequency
#' @importFrom stats reshape
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of stacked bar graphs of amino acid proportions
percentAA <- function(input.data, 
                      chain = "TRB", 
                      group.by = NULL, 
                      order.by = NULL,
                      aa.length = 20,
                      exportTable = FALSE, 
                      palette = "inferno",
                      ...)  {
  
  sco <- .is.seurat.or.se.object(input.data)
  input.data <- .dataWrangle(input.data, group.by, "CTaa", chain)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  # Calculating proportion 
  lapply(seq_along(input.data), function(x) {
    strings <- .processStrings(input.data[[x]][["CTaa"]], aa.length)
    freqs <- calculateFrequency(strings)    
    freqs_df <- as.data.frame(freqs)
    freqs_df$AminoAcid <- rownames(freqs_df)
    long_df <- reshape(
      freqs_df,
      varying = list(colnames(freqs)), 
      v.names = "Frequency",           
      timevar = "Position",            
      times = colnames(freqs),         
      direction = "long",              
      idvar = "AminoAcid"              
    )
    rownames(long_df) <- NULL # Remove the messy row names
    long_df$Position <- as.integer(gsub("Pos.", "", long_df$Position))
    long_df$AminoAcid[long_df$AminoAcid == "."] <- NA
    long_df$group <- names(input.data)[x]
    long_df
  }) -> res.list
  
  mat_melt <- do.call(rbind, res.list)
  if(!is.null(order.by)) {
    mat_melt <- .orderingFunction(vector = order.by,
                                  group.by = "group", 
                                  mat_melt)
  }
  if (exportTable == TRUE) { 
    return(mat_melt) 
  }
  
  # Plotting the result
  plot <- ggplot(mat_melt, aes(x=as.factor(Position), y = .data[["Frequency"]], fill=.data[["AminoAcid"]])) +
    geom_bar(stat = "identity", position="fill", lwd= 0.25, color = "black") +
    scale_fill_manual(name = "Amino Acid", 
                      values = rev(.colorizer(palette,20))) +
    xlab("Amino Acid Residues") +
    ylab("Relative Percent") +
    .themeRepertoire(...) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if(length(res.list) > 1) {
    plot <- plot + facet_grid(group~.)
  }
  return(plot)
}    
    
