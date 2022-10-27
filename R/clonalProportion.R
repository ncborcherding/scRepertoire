#' Examining the clonal space occupied by specific clonotypes
#'
#' This function calculates the relative clonal space occupied by the 
#' clonotypes. The grouping of these clonotypes is based on the parameter 
#' split, at default, split will group the clonotypes into bins of 1:10, 
#' 11:100, 101:1001, etc. To adjust the clonotypes selected, change the 
#' numbers in the variable split. If a matrix output for the data is
#' preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param split The cutpoints for the specific clonotypes.
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
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#'
#' @export
#' @return ggplot of the space occupied by the specific rank of clonotypes
clonalProportion <- function(df,
                             split = c(10, 100, 1000, 10000, 30000, 100000), 
                             cloneCall = "strict", 
                             chain = "both", 
                             group.by = NULL,
                             split.by = NULL,
                             exportTable = FALSE) {
    Con.df <- NULL
    df <- list.input.return(df, split.by = split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkList(df)
    df <- checkBlanks(df, cloneCall)
    if(!is.null(group.by)) {
      df <- groupList(df, group.by)
    }
    mat <- matrix(0, length(df), length(split), dimnames = list(names(df), 
            paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))
    if (chain != "both") {
      for (x in seq_along(df)) {
        df[[x]] <- off.the.chain(df[[x]], chain, cloneCall)
      }
    }
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, split[-length(split)] + 1)
    for (i in seq_along(split)) {
        mat[,i] <- vapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])), 
                            FUN.VALUE = numeric(1))
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", 
                        values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()
    return(plot)



}
