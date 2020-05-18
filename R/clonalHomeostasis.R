#' Examining the clonal homeostasis
#'
#' This function calculates the space occupied by specific clonotype proportions. The grouping of these clonotypes is
#' based on the parameter cloneTypes, at default, cloneTypes will group the clonotypes into bins of Rare = 0 to 0.0001,
#' Small = 0.0001 to 0.001, etc. To adjust the proportions, change the number or labeling of the cloneTypes paramter. If a
#' matrix output for the data is preferred, set exportTable = TRUE.
#'
#' @param df The product of CombineContig() or expression2List()
#' @param cloneTypes The cutpoints of the proportions.
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) or CDR3 amino acid (aa), or
#' CDR3 gene+nucleotide (gene+nt).
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @export
clonalHomeostasis <- function(df,
                              cloneTypes = c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1),
                              cloneCall = c("gene", "nt", "aa", "gene+nt"),
                              exportTable = F) {
    cloneTypes <- c(None = 0, cloneTypes)

    cloneCall <- theCall(cloneCall)
    df <- if(is(df)[1] != "list") list(df) else df

    mat <- matrix(0, length(df), length(cloneTypes) - 1, dimnames = list(names(df), names(cloneTypes)[-1]))
    df <- lapply(df, '[[', cloneCall)
    for (i in seq_along(df)) {
        df[[i]] <- na.omit(df[[i]])
    }

    fun <- function(x) {
        table(x)/length(x)
    }

    df <- lapply(df, fun)

    for (i in 2:length(cloneTypes)) {
        mat[,i-1] <- sapply(df, function (x) sum(x[x > cloneTypes[i-1] & x <= cloneTypes[i]]))
        colnames(mat)[i-1] <- paste0(names(cloneTypes[i]), ' (', cloneTypes[i-1], ' < X <= ', cloneTypes[i], ')')
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonotype Group", values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Relative Abundance") +
        theme_classic()

    return(plot)

}
