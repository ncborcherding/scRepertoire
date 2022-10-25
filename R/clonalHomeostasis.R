#' Examining the clonal homeostasis
#'
#' This function calculates the space occupied by clonotype proportions. 
#' The grouping of these clonotypes is based on the parameter cloneTypes, 
#' at default, cloneTypes will group the clonotypes into bins of Rare = 0 
#' to 0.0001, Small = 0.0001 to 0.001, etc. To adjust the proportions, 
#' change the number or labeling of the cloneTypes paramter. If a matrix 
#' output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalHomeostasis(combined, cloneCall = "gene")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneTypes The cutpoints of the proportions.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#' @export
#' @return ggplot of the space occupied by the specific proportion of clonotypes
clonalHomeostasis <- function(df, 
                              cloneTypes = c(Rare = .0001, Small = .001, 
                              Medium = .01, Large = .1, Hyperexpanded = 1),
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

    mat <- matrix(0, length(df), length(cloneTypes) - 1, 
                dimnames = list(names(df), 
                names(cloneTypes)[-1]))
    if (chain != "both") {
      for (x in seq_along(df)) {
        df[[x]] <- off.the.chain(df[[x]], chain, cloneCall)
      }
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
        scale_fill_manual(name = "Clonotype Group", 
                    values = rev(colorblind_vector(col))) +
        xlab("Samples") +
        ylab("Relative Abundance") +
        theme_classic()
    return(plot)
}
