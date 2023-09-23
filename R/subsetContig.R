
#' Subset the product of combineTCR() combineBCR()
#'
#' This function allows for the subsetting of the product of 
#' \code{\link{combineTCR}} or \code{\link{combineBCR}} 
#' by the name of the individual list element. 
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' subset <- subsetContig(combined, name = "sample", variables = c("P17B"))
#'
#' @param df The product of \code{\link{combineTCR}} or \code{\link{combineBCR}}.
#' @param name The column header you'd like to use to subset.
#' @param variables The values to subset by, must be in the names(df).

#' @export
#' @return list of contigs that have been filtered for the name parameter
subsetContig <- function(df, name, variables = NULL) {
  names2 <- NULL
  df2 <- list()
  for (i in seq_along(df)) {
    if (df[[i]][1,name] %in% variables) {
      df2 <- append(df2, list(df[[i]]))
      n2 <- names(df)[i] 
      names2 <- c(names2, n2)
    }
    else {
      next()
    }
  }
  names(df2) <- names2
  return(df2)
}