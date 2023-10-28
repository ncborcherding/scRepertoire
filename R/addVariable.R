#' Adding variables after the combination of clonotypes.
#'
#' This function adds variables to the product of \code{\link{combineTCR}}, 
#' or \code{\link{combineBCR}} to be used in later visualizations. 
#' For each element, the function will add a column (labeled by name) 
#' with the variable. The length of the variable paramater needs to match
#' the length of the combined object.
#'
#' @examples
#' combined <- combineTCR(contig_list, 
#'                        samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' combined <- addVariable(combined, 
#'                         variable.name = "Type", 
#'                         variables = rep(c("B", "L"), 4))
#'
#' @param input.data The product of \code{\link{combineTCR}} or \code{\link{combineBCR}}.
#' @param variable.name The column header to add.
#' @param variables The exact values to add to each element of the list.
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return input.data list with the variable column added to each element.
addVariable <- function(input.data, 
                        variable.name = NULL, 
                        variables =  NULL) {
  if (length(input.data) != length(variables)) {
    stop("Make sure the variables match the length of the contig list")
  }
  for (i in seq_along(input.data)) {
    input.data[[i]][,variable.name] <- variables[i]
  }
  return(input.data)
}