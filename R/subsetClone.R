#' Subset The Product of combineTCR() or combineBCR()
#'
#' This function allows for the subsetting of the product of 
#' [combineTCR()] or [combineBCR()] 
#' by the name of the individual list element. 
#'
#' @examples
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' subset <- subsetClones(combined, name = "sample", variables = c("P17B"))
#'
#' @param input.data The product of [combineTCR()] or 
#' [combineBCR()].
#' @param name The column header/name to use for subsetting.
#' @param variables The values to subset by, must be in the names(input.data).

#' @export
#' @concept Loading_and_Processing_Contigs
#' @return list of contigs that have been filtered for the name parameter
subsetClones <- function(input.data, 
                         name, 
                         variables = NULL) {
  names2 <- NULL
  input.data2 <- list()
  for (i in seq_along(input.data)) {
    if (input.data[[i]][1,name] %in% variables) {
      input.data2 <- append(input.data2, list(input.data[[i]]))
      n2 <- names(input.data)[i] 
      names2 <- c(names2, n2)
    }
    else {
      next()
    }
  }
  names(input.data2) <- names2
  return(input.data2)
}
