#' Generate data frame to be used with circlize R package to visualize
#' clonotypes as a chord diagram. 
#' 
#' This function will take the meta data from the product of 
#' \code{\link{combineExpression}} and generate a relational data frame to 
#' be used for a chord diagram. Each cord will represent the number of 
#' clonotype unique and shared across the multiple \strong{group.by} variable. 
#' If using the downstream circlize R package, please read and cite the
#' following \href{https://pubmed.ncbi.nlm.nih.gov/24930139/}{manuscript}.
#' If looking for more advance ways for circular visualizations, there
#' is a great \href{https://jokergoo.github.io/circlize_book/book/}{cookbook}
#' for the circlize package.
#' 
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined, 
#'                                    scRep_example)
#' 
#' #Getting data frame output for Circlize
#' circles <- getCirclize(scRep_example, 
#'                        group.by = "seurat_clusters")
#' 
#' 
#' @param sc.data The single-cell object after \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param proportion Boolean will calculate relationship unique 
#' clonotypes (proportion = FALSE) or a ratio of the group.by 
#' variable (proportion = TRUE)
#' 
#' @importFrom reshape2 dcast
#' @export
#' @concept SC_Functions
#' @return data frame of shared clonotypes between groups
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(sc.data, 
                        cloneCall = "strict", 
                        group.by = NULL, 
                        proportion = FALSE) {
  meta <- .grabMeta(sc.data)
  cloneCall <- .theCall(meta, cloneCall)
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  #Quantifying clones across group.by variable
  dat <- meta[, c(cloneCall, group.by)]
  dat <- dat[!is.na(dat[,cloneCall]),]
  mat <- suppressMessages(dcast(dat, dat[,cloneCall] ~ dat[,group.by]))
  mat <- mat[apply(mat[,-1], 1, function(x) !all(x==0)),]
  mat2<- mat[-1]
  mat2[mat2 >= 1] <- 1
  total <- nrow(mat)
  
  list <- list()
  for (i in seq_len(nrow(mat2))) {
    list[[i]] <- which(mat2[i,] > 0)
  }
  
  #Making a pairwise matrix that sums the shared clonotypes
  matrix_out <- matrix(ncol = ncol(mat2), nrow = ncol(mat2), 0)
  for (j in seq_along(list)) {
    matrix_out[list[[j]],list[[j]]] <- matrix_out[list[[j]],list[[j]]] + 1
    if (length(list[[j]]) > 1) {
      #length <- length(list[[j]])
      diag(matrix_out[list[[j]],list[[j]]]) <-  diag(matrix_out[list[[j]],list[[j]]]) - 1
    }
  }
  matrix_out[lower.tri(matrix_out)] <- NA
  
  colnames(matrix_out) <- colnames(mat2)
  rownames(matrix_out) <- colnames(mat2)
  
  #Converting matrix to data frame to work with criclize
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
