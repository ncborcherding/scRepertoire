
  
#' @param df The product of combineTCR(), combineBCR(), 
#' expression2List(), or combineExpression().
#' @param chain "TRA", "TRB", "TRG", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param order Maintain the order of the list when plotting
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param exportTable Returns the data frame used for forming the graph
#' @import ggplot2
#' @export
#' @return ggplot of the total or relative unique clonotypes
percentAA <- function(df, 
                        chain = "TRB", 
                        group.by = NULL, 
                        split.by = NULL,
                        exportTable = FALSE) 
  df <- list.input.return(df, split.by)
  df <- checkBlanks(df, "CTaa")
  for(i in seq_along(df)) {
      df[[i]] <- off.the.chain(df[[i]], chain, "CTaa")
  }
  
  for (i in seq_along(df)) {
    strings <- df[[i]][,"CTaa"]
    strings <- do.call(c,str_split(strings, ";"))
    strings <- strings[strings != "NA"]
    aa.length <- max(nchar(strings))
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
    ######
    #Need to come up with a summary function
  }
  #Need to come up with a visualization
    
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
