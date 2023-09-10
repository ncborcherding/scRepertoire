#' Adding variables after the combination of contigs.
#'
#' This function adds variables to the product of combineTCR() combineBCR() or 
#' expression2List() to be used in later visualizations. For each element, 
#' the function will add a column (labled by name) with the variable. 
#' The length of the variable paramater needs to match the length of 
#' the combined object.
#'
#' @examples
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' combined <- addVariable(combined, 
#'                        name = "Type", 
#'                        variables = rep(c("B", "L"), 4))
#'
#' @param df The product of combineTCR() combineBCR() or expression2List().
#' @param name The column header to add.
#' @param variables The exact values to add to each element of the list.
#' @export
#' @return list of contigs with a new column (name).
addVariable <- function(df, name = NULL, variables =  NULL) {
  if (length(df) != length(variables)) {
    stop("Make sure the variables match the length of the contig list")
  }
  for (i in seq_along(df)) {
    df[[i]][,name] <- variables[i]
  }
  return(df)
}