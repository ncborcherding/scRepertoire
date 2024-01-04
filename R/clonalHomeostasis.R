#' Examining the clonal homeostasis
#'
<<<<<<< Updated upstream
#' This function calculates the space occupied by clonotype proportions. 
#' The grouping of these clonotypes is based on the parameter cloneTypes, 
#' at default, cloneTypes will group the clonotypes into bins of Rare = 0 
=======
#' This function calculates the space occupied by clone proportions. 
#' The grouping of these clones is based on the parameter cloneSize, 
#' at default, cloneSize will group the clones into bins of Rare = 0 
>>>>>>> Stashed changes
#' to 0.0001, Small = 0.0001 to 0.001, etc. To adjust the proportions, 
#' change the number or labeling of the cloneTypes paramter. If a matrix 
#' output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' clonalHomeostasis(combined, cloneCall = "gene")
#'
<<<<<<< Updated upstream
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneTypes The cutpoints of the proportions.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
=======
#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}.
#' @param cloneSize The cut points of the proportions.
#' @param cloneCall How to call the clone - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa),
#' VDJC gene + CDR3 nucleotide (strict) or a custom variable in the data. 
>>>>>>> Stashed changes
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
<<<<<<< Updated upstream
=======
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}.
>>>>>>> Stashed changes
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#' @export
<<<<<<< Updated upstream
#' @return ggplot of the space occupied by the specific proportion of clonotypes
clonalHomeostasis <- function(df, 
                              cloneTypes = c(Rare = .0001, Small = .001, 
                              Medium = .01, Large = .1, Hyperexpanded = 1),
=======
#' @concept Visualizing_Clones
#' @return ggplot of the space occupied by the specific proportion of clones
clonalHomeostasis <- function(input.data, 
                              cloneSize = c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1),
>>>>>>> Stashed changes
                              cloneCall = "strict", 
                              chain = "both", 
                              group.by = NULL,
                              split.by = NULL,
                              exportTable = FALSE) {
    cloneTypes <- c(None = 0, cloneTypes)
    df <- list.input.return(df, split.by = split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkList(df)
    df <- checkBlanks(df, cloneCall)
    if(!is.null(group.by)) {
      df <- groupList(df, group.by)
    }

<<<<<<< Updated upstream
    mat <- matrix(0, length(df), length(cloneTypes) - 1, 
                dimnames = list(names(df), 
                names(cloneTypes)[-1]))
    if (chain != "both") {
      for (x in seq_along(df)) {
        df[[x]] <- off.the.chain(df[[x]], chain, cloneCall)
      }
=======
    #Assigning the clonal grouping
    input.data <- lapply(input.data, '[[', cloneCall)
    input.data <- lapply(input.data, na.omit)
    fun <- function(x) { table(x)/length(x) }
    input.data <- lapply(input.data, fun)
    for (i in 2:length(cloneSize)) {
        mat[,i-1] <- vapply(input.data, function (x) 
                              sum(x[x > cloneSize[i-1] & x <= cloneSize[i]]), FUN.VALUE = numeric(1))
        colnames(mat)[i-1] <- paste0(names(cloneSize[i]), ' (', 
                                    cloneSize[i-1], ' < X <= ', 
                                    cloneSize[i], ')') }
    if (exportTable) { 
      return(mat) 
>>>>>>> Stashed changes
    }
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    fun <- function(x) { table(x)/length(x) }
    df <- lapply(df, fun)
    for (i in 2:length(cloneTypes)) {
        mat[,i-1] <- vapply(df, function (x) sum(x[x > cloneTypes[i-1] & x <= 
                            cloneTypes[i]]), FUN.VALUE = numeric(1))
        colnames(mat)[i-1] <- paste0(names(cloneTypes[i]), ' (', 
                                    cloneTypes[i-1], ' < X <= ', 
                                    cloneTypes[i], ')') }
    if (exportTable == TRUE) { return(mat) }
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
<<<<<<< Updated upstream
        scale_fill_manual(name = "Clonotype Group", 
                    values = rev(colorblind_vector(col))) +
=======
        scale_fill_manual(name = "Clonal Group", 
                    values = .colorizer(palette,col)) +
>>>>>>> Stashed changes
        xlab("Samples") +
        ylab("Relative Abundance") +
        theme_classic()
    return(plot)
}
