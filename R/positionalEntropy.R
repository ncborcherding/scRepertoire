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
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}.
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL".
#' @param group.by The variable to use for grouping.
#' @param aa.length The maximum length of the CDR3 amino acid sequence. 
#' @param method The method to calculate the entropy/diversity - 
#' "shannon", "inv.simpson", "norm.entropy".
#' @param n.boots number of bootstraps to down sample in order to 
#' get mean diversity.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom stringr str_split
#' @export
#' @concept Summarize_Repertoire
#' @return ggplot of line graph of diversity by position
positionalEntropy <- function(input.data, 
                              chain = "TRB", 
                              group.by = NULL, 
                              aa.length = 20,
                              method = "shannon",
                              n.boots = 20,
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
  
  #Selecting Diversit Function
  diversityFunc <- switch(method,
                          "norm.entropy" = .shannon,
                          "inv.simpson" = .invsimpson,
                          "shannon" = .normentropy,
                          stop("Invalid method provided"))
  
  min <- .short.check(input.data, cloneCall)
  
  lapply(input.data, function(x) {
      lapply(seq_len(n.boots), function(y) {
       strings <- x[,cloneCall]
       strings <- do.call(c,str_split(strings, ";"))
       strings <- strings[strings != "NA"]
       strings <- na.omit(strings)
       strings <- strings[nchar(strings) < aa.length]
       strings <- strings[sample(seq_len(length(strings)), min)]
       strings <- .padded_strings(strings, aa.length)
       strings <- do.call(rbind, strings)
       aa.output <- apply(strings, 2, function(z) {
         summary <- as.data.frame(table(z, useNA = "always"))
       })
       res <- suppressWarnings(Reduce(function(...) merge(..., all = TRUE, by="z"), aa.output))
       colnames(res) <- c("AA", paste0("pos.", seq_len(aa.length)))
       res[seq_len(20),][is.na(res[seq_len(20),])] <- 0
       diversity <- sapply(res[,2:ncol(res)], diversityFunc)
       diversity[is.nan(diversity)] <- 0
       diversity
    }) -> diversity.calculations
    diversity.calculations <- do.call(rbind, diversity.calculations)
    diversity.means <- colMeans(diversity.calculations)
    diversity.means
    }) -> positional.diversity
    
    mat <- do.call(rbind, positional.diversity)
    mat_melt <- suppressMessages(melt(mat))
    
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

  


