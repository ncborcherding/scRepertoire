
#' Subset the product of combineTCR() combineBCR() or expression2List()
#'
#' This function allows for the subsetting of the product of combineTCR() 
#' combineBCR() or expression2List() by the name of the individual list 
#' element. In general the names of are samples + _ + ID, allowing 
#' for users to subset the product of combineTCR(), combineBCR(), 
#' or expression2List() across a string or individual name.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' subset <- subsetContig(combined, name = "sample", variables = c("PX"))
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List().
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