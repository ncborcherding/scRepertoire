#' Examining the clonal space occupied by specific clonotypes
#'
#' @description
#' Calculate the relative clonal space occupied by the clonotypes. To adjust the clonotypes selected, change the numbers in the variable split.
#'
#' @param df The product of CombineContig() or the seurat object after combineSeurat()
#' @param split The cutpoints for the specific clonotypes
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#'
#' @export
clonalProportion <- function(df,
                             split = c(10, 100, 1000, 10000, 30000, 100000),
                             cloneCall = c("gene", "nt", "aa", "gene+nt"),
                             exportTable = F) {
    Con.df <- NULL
    cloneCall <- theCall(cloneCall)

    if (inherits(x=df, what ="Seurat")) {
        meta <- data.frame(df@meta.data, df@active.ident)
        colnames(meta)[length(meta)] <- "cluster"
        unique <- str_sort(as.character(unique(meta$cluster)), numeric = TRUE)
        df <- NULL
        for (i in seq_along(unique)) {
            subset <- subset(meta, meta[,"cluster"] == unique[i])
            df[[i]] <- subset
        }
        names(df) <- unique
    }
    df <- if(class(df) != "list") list(df) else df

    mat <- matrix(0, length(df), length(split), dimnames = list(names(df), paste0('[', c(1, split[-length(split)] + 1), ':', split, ']')))
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- na.omit(df[[i]])
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, split[-length(split)] + 1)
    for (i in 1:length(split)) {
        mat[,i] <- sapply(df, function (x) sum(na.omit(x[cut[i]:split[i]])))
    }
    if (exportTable == T) {
        clonalProportion_output <<- mat
    }
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()
    return(plot)



}
