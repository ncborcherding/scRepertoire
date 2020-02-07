#' Examining the clonal homeostasis
#'
#' @param df The product of CombineContig() or the seurat object after combineSeurat()
#' @param cloneTypes The cutpoints of the proportions
#' @param call How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#'
#' @export
clonalHomeostasis <- function(df,
                              cloneTypes = c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1),
                              call = c("gene", "nt", "aa", "gene+nt")) {
    require(ggplot2)
    cloneTypes <- c(None = 0, cloneTypes)
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else if (call == "aa") {
        call <- "CTaa"
    } else if (call == "gene+nt") {
        call <- "CTstrict"
    } else {
        stop("Are you sure you made the right call? ", .call = F)
    }
    if (class(df)[1] == "Seurat") {
        meta <- data.frame(df@meta.data, df@active.ident)
        colnames(meta)[length(meta)] <- "cluster"
        unique <- stringr::str_sort(as.character(unique(meta$cluster)), numeric = TRUE)
        df <- NULL
        for (i in seq_along(unique)) {
            subset <- subset(meta, meta[,"cluster"] == unique[i])
            df[[i]] <- subset
        }
        names(df) <- unique
    }
    df <- if(class(df) != "list") list(df) else df

    mat <- matrix(0, length(df), length(cloneTypes) - 1, dimnames = list(names(df), names(cloneTypes)[-1]))
    df <- lapply(df, '[[', call)
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

    mat_melt <- reshape2::melt(mat)
    col <- length(unique(mat_melt$Var2))
    ggplot2::ggplot(mat_melt, aes(x=Var1, y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonotype Group", values = colorblind_vector(col)) +
        xlab("Samples") +
        ylab("Relative Abundance") +
        theme_classic()

}
