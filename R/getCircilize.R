#' Generate data frame to be used with circlize R package to visualize
#' clonotypes as a chord diagram. 
#' 
#' This function will take the meta data from the product of 
#' combineExpression()and generate a relational data frame to 
#' be used for a chord diagram. Each cord will represent the number of 
#' clonotype unqiue and shared across the multiple group.by variable.
#' 
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' #Getting data frame output for Circilize
#' circles <- getCirclize(scRep_example, group.by = "seurat_clusters")
#' 
#' 
#' @param sc The single-cell object after \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param proportion Binary will calculate relationship unique 
#' clonotypes (proportion = FALSE) or a ratio of the group.by 
#' variable (proportion = TRUE)
#' 
#' @importFrom reshape2 dcast
#' @export
#' @return data frame of shared clonotypes between groups
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(sc, cloneCall = "strict", 
                        group.by = NULL, proportion = FALSE) {
  meta <- grabMeta(sc)
  cloneCall <- theCall(cloneCall)
  test <- meta[, c(cloneCall, group.by)]
  test <- test[!is.na(test[,cloneCall]),]
  dTest <- suppressMessages(dcast(test, test[,cloneCall] ~ test[,group.by]))
  dTest <- dTest[apply(dTest[,-1], 1, function(x) !all(x==0)),]
  dTest2 <- dTest[-1]
  dTest2[dTest2 >= 1] <- 1
  total <- nrow(dTest)
  
  list <- list()
  for (i in seq_len(nrow(dTest2))) {
    list[[i]] <- which(dTest2[i,] > 0)
  }
  matrix_out <- matrix(ncol = ncol(dTest2), nrow = ncol(dTest2), 0)
  for (j in seq_along(list)) {
    matrix_out[list[[j]],list[[j]]] <- matrix_out[list[[j]],list[[j]]] + 1
    if (length(list[[j]]) > 1) {
      #length <- length(list[[j]])
      diag(matrix_out[list[[j]],list[[j]]]) <-  diag(matrix_out[list[[j]],list[[j]]]) - 1
    }
  }
  
  matrix_out[lower.tri(matrix_out)] <- NA
  
  colnames(matrix_out) <- colnames(dTest2)
  rownames(matrix_out) <- colnames(dTest2)
  
  output <- data.frame(from = rep(rownames(matrix_out), 
                                  times = ncol(matrix_out)),
                       to = rep(colnames(matrix_out), each = nrow(matrix_out)),
                       value = as.vector(matrix_out),
                       stringsAsFactors = FALSE)
  output <- na.omit(output)
  
  if (proportion == TRUE) {
    output$value <- output$value/total
  } 
  return(output)
}