#' Examining the relative amino acid composition by position
#'
#' This function the proportion of amino acids along the residues 
#' of the cdr3 amino acid sequence.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentAA(combined, chain = "TRB", aa.length = 20)
  
#' @param df The product of combineTCR(), combineBCR(), 
#' expression2List(), or combineExpression().
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param aa.length The maximum length of the cdr3 amino acid sequence 
#' @param exportTable Returns the data frame used for forming the graph
#' @import ggplot2
#' @export
#' @return ggplot of stacked bar graphs of amino acid proportions
percentAA <- function(df, 
                        chain = "TRB", 
                        group.by = NULL, 
                        split.by = NULL,
                        aa.length = 30,
                        exportTable = FALSE) 
  df <- list.input.return(df, split.by)
  df <- checkBlanks(df, "CTaa")
  for(i in seq_along(df)) {
      df[[i]] <- off.the.chain(df[[i]], chain, "CTaa")
  }
  res.list <- list()
  for (i in seq_along(df)) {
    strings <- df[[i]][,"CTaa"]
    strings <- do.call(c,str_split(strings, ";"))
    strings <- strings[strings != "NA"]
    strings <- strings[nchar(strings) < aa.length]
    strings <- .padded_strings(strings, aa.length)
    strings <- do.call(rbind, strings)
    
    #Summarizing the % of each position
    aa.output <- apply(strings, 2, function(x) {
      summary <- as.data.frame(prop.table(table(x, useNA = "always")))
    })
    
    #Forming a matrix of % across each position and formatting
    res <- suppressWarnings(Reduce(function(...) merge(..., all = TRUE, by="x"), aa.output))
    colnames(res) <- c("AA", paste0("pos.", seq_len(aa.length)))
    res[1:20,][is.na(res[1:20,])] <- 0
    melt.res <- suppressMessages(melt(res))
    melt.res$group <- names(df)[i]
    res.list[[i]] <- melt.res
  }
  mat_melt <- do.call(rbind, res.list)
  plot <- ggplot(mat_melt, aes(x=as.factor(variable), y = value, fill=AA)) +
    geom_bar(stat = "identity", position="fill", lwd= 0.25, color = "black") +
    scale_fill_manual(name = "Amino Acid", 
                      values = rev(colorblind_vector(21))) +
    xlab("Amino Acid Residues") +
    ylab("Relative Percent") +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if(length(res.list) > 1) {
    plot <- plot + facet_grid(group~.)
  }
  if (exportTable == TRUE) { return(mat_melt) }
  return(plot)
}    
    
.padded_strings <- function(strings, max.length) {
      max_length <- max.length
      
      x <- lapply(strings, function(str) {
        str_len <- nchar(str)
        str <- strsplit(str, split = "")[[1]]
        if (str_len < max_length) {
          c(str, rep(NA, max_length - str_len))
        } else {
          str
        }
    })
  }
