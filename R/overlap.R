#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the caclulation and visualizations of the overlap coefficient or morisita index
#' for clonotypes using the product of combineContig() or expression2list(). The overlap coefficient is calculated
#' using the intersection of clonotypes divided by the length of the smallest componenet. Morisita index is estimating
#' the dispersion of a population, more information can be found [here](https://en.wikipedia.org/wiki/Morisita%27s_overlap_index).
#' If a matrix output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' clonalOverlap(combined, cloneCall = "gene", method = "overlap")
#'
#' @param df The product of CombineContig() or expression2List()
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) or CDR3 amino acid (aa), or
#' CDR3 gene+nucleotide (gene+nt).
#' @param method The method to calculate the overlap, either the overlap coefficient or morisita index
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list
clonalOverlap <- function(df,
                    cloneCall = c("gene", "nt", "aa", "gene+nt"),
                    method = c("overlap", "morisita"),
                    exportTable = FALSE){

    cloneCall <- theCall(cloneCall)

    df <- df[order(names(df))]
    values <- str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]

    num_samples <- length(df[])
    names_samples <- names(df)
    length <- 1:num_samples
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples
    if (method == "overlap") {

        for (i in seq_along(length)){
          df.i <- df[[i]]
          df.i <- df.i[,c("barcode",cloneCall)]
          df.i_unique <- df.i[!duplicated(df.i$barcode),]

          for (j in seq_along(length)){
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

        for (i in seq_along(length)){
            df.i <- df[[i]]
            df.i <- data.frame(table(df.i[,cloneCall]))
            colnames(df.i) <- c(cloneCall, 'Count')
            df.i[,2] <- as.numeric(df.i[,2])

            for (j in seq_along(length)){
                if (i >= j){
                    next
                }
                else {
                    df.j <- df[[j]]
                    df.j <- data.frame(table(df.j[,cloneCall]))
                    colnames(df.j) <- c(cloneCall, 'Count')
                    df.j[,2] <- as.numeric(df.j[,2])

                    merged <- merge(df.i, df.j, by = cloneCall, all = TRUE)
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
