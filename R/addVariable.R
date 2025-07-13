#' Adding Variables After combineTCR() or combineBCR()
#'
#' This function adds variables to the product of [combineTCR()], 
#' or [combineBCR()] to be used in later visualizations. 
#' For each element, the function will add a column (labeled by 
#' **variable.name**) with the variable. The length of the 
#' **variables** parameter needs to match the length of the 
#' combined object.
#'
#' @examples
#' combined <- combineTCR(contig_list, 
#'                        samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' combined <- addVariable(combined, 
#'                         variable.name = "Type", 
#'                         variables = rep(c("B", "L"), 4))
#'
#' @param input.data The product of [combineTCR()] or [combineBCR()].
#' @param variable.name A character string that defines the new variable to add.
#' @param variables A character vector defining the desired column value for 
#' each list element. Must match the length of the product of [combineTCR()] or 
#' [combineBCR()].
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
