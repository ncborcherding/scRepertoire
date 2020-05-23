#' Examining the clonal space occupied by specific clonotypes
#'
#' This function calculates the relative clonal space occupied by the clonotypes. 
#' The grouping of these clonotypes is based on the parameter split, at default, 
#' split will group the clonotypes into bins of 1:10, 11:100, 101:1001, etc. 
#' To adjust the clonotypes selected, change the numbers in the variable split. 
#' If a matrix output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), rep(c("P", "T"), 3), cells ="T-AB")
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param df The product of CombineContig() or expression2List()
#' @param split The cutpoints for the specific clonotypes.
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) 
#' or CDR3 amino acid (aa), or CDR3 gene+nucleotide (gene+nt).
#' @param exportTable Exports a table of the data into the global environment in 
#' addition to the visualization
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#'
#' @export
#' @return ggplot of the space occupied by the specific rank of clonotypes
clonalProportion <- function(df,
                                split = c(10, 100, 1000, 10000, 30000, 100000),
                                cloneCall = c("gene", "nt", "aa", "gene+nt"),
                                exportTable = FALSE) {
    Con.df <- NULL
    cloneCall <- theCall(cloneCall)
    df <- if(is(df)[1] != "list") list(df) else df

    mat <- matrix(0, length(df), length(split), dimnames = list(names(df), 
                                                                paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- na.omit(df[[i]])
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, split[-length(split)] + 1)
    for (i in seq_along(split)) {
        mat[,i] <- vapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])), FUN.VALUE = numeric(1))
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()
    return(plot)



}
