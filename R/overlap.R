#' Examining the clonal overlap between samples
#'
#' @description
#' Use overlap coefficient or morisita index to exam the overlapping clonotypes across samples. If a matrix output for the data is preferred, set exportTable = T.
#'
#' @param df The product of CombineContig() or the seurat object after combineSeurat()
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param method The method to calculate the overlap, either the overlap coefficient or morisita index
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @export
clonalOverlap <- function(df,
                    cloneCall = c("gene", "nt", "aa", "gene+nt"),
                    method = c("overlap", "morisita"),
                    exportTable = F){

    cloneCall <- theCall(cloneCall)

    if (is(df)[1] == "Seurat") {
        meta <- data.frame(df@meta.data, df@active.ident)
        colnames(meta)[ncol(meta)] <- "cluster"
        unique <- str_sort(as.character(unique(meta$cluster)), numeric = TRUE)
        meta$barcode <- rownames(meta)
        df <- NULL
        for (i in seq_along(unique)) {
            subset <- subset(meta, meta[,"cluster"] == unique[i])
            df[[i]] <- subset
        }
        names(df) <- unique
    }
    else if (is(df)[1] != "Seurat") {
        df <- df[order(names(df))]
    }

    num_samples <- length(df[])
    names_samples <- names(df)
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples
    if (method == "overlap") {

        for (i in 1:num_samples){
          df.i <- df[[i]]
          df.i <- df.i[,c("barcode",cloneCall)]
          df.i_unique <- df.i[!duplicated(df.i$barcode),]

          for (j in 1:num_samples){
            if (i >= j){
              next
            }
            else {
              df.j <- df[[j]]
              df.j <- df.j[,c("barcode",cloneCall)]
              df.j_unique <- df.j[!duplicated(df.j$barcode),]
              overlap <- length(intersect(df.i_unique[,cloneCall], df.j_unique[,cloneCall]))
              coef_matrix[i,j] <- overlap/min(length(df.i_unique[,cloneCall]), length(df.j_unique[,cloneCall]))
            }
          }
        }
    }
    else if (method == "morisita") {

        for (i in 1:num_samples){
            df.i <- df[[i]]
            df.i <- data.frame(table(df.i[,cloneCall]))
            colnames(df.i) <- c(cloneCall, 'Count')
            df.i[,2] <- as.numeric(df.i[,2])

            for (j in 1:num_samples){
                if (i >= j){
                    next
                }
                else {
                    df.j <- df[[j]]
                    df.j <- data.frame(table(df.j[,cloneCall]))
                    colnames(df.j) <- c(cloneCall, 'Count')
                    df.j[,2] <- as.numeric(df.j[,2])

                    merged <- merge(df.i, df.j, by = cloneCall, all = T)
                    merged[is.na(merged)] <- 0
                    sum.df.i <- sum(df.i[,2])
                    sum.df.j <- sum(df.j[,2])
                    coef.i.j <- 2 * sum(merged[,2] * merged[,3] / sum.df.j) / sum.df.j/
                        ((sum((df.i[,2] / sum.df.i)^2) + sum((df.j[,2] / sum.df.j)^2)))

                    coef_matrix[i,j] <- coef.i.j
                }
            }
        }

    }
    coef_matrix$names <- rownames(coef_matrix)
    if (exportTable == TRUE) {
        return(coef_matrix)
    }
    coef_matrix <- suppressMessages(melt(coef_matrix))
    coef_matrix <- coef_matrix[,-1]
    col <- colorblind_vector(7)
    values <- str_sort(as.character(unique(coef_matrix$names)), numeric = TRUE)
    values2 <- quiet(dput(values))
    coef_matrix$variable <- factor(coef_matrix$variable, levels = values2)
    coef_matrix$names <- factor(coef_matrix$names, levels = values2)
    plot <- ggplot(coef_matrix, aes(x=names, y=variable, fill=value)) +
        geom_tile() +
        geom_text(aes(label = round(value, digits = 3))) +
        scale_fill_gradient2(high = col[1], mid = col[4], midpoint = ((range(na.omit(coef_matrix$value)))/2)[2], low=col[7], na.value = "white") +
        labs(fill = method) +
        theme_classic() +
        theme(axis.title = element_blank())
    return(plot)
}
