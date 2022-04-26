#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of the 
#' overlap coefficient, morisita, or jaccard index for clonotypes 
#' using the product of combineTCR(), combineBCR() or expression2list(). 
#' The overlap coefficient is calculated using the intersection of clonotypes 
#' divided by the length of the smallest component. 
#' If a matrix output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' clonalOverlap(combined, cloneCall = "gene", method = "overlap")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the overlap, either the "overlap" 
#' coefficient, "morisita", "jaccard" indices, or "raw" for the base numbers.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list
clonalOverlap <- function(df, cloneCall = "strict", 
                                method = c("overlap", "morisita", "jaccard", "raw"), 
                                chain = "both", 
                                split.by = NULL, 
                                exportTable = FALSE){
    df <- list.input.return(df, split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    df <- df[order(names(df))]
    values <- str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]
    num_samples <- length(df[])
    names_samples <- names(df)
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples
    length <- seq_len(num_samples)
    if (chain != "both") {
      for (i in seq_along(df)) {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
    }
    if (method == "overlap") {
        coef_matrix <- overlapIndex(df, length, cloneCall, coef_matrix)
    } else if (method == "morisita") {
        coef_matrix <- morisitaIndex(df, length, cloneCall, coef_matrix)
    } else if (method == "jaccard") {
        coef_matrix <- jaccardIndex(df, length, cloneCall, coef_matrix)
    } else if (method == "raw") {
        coef_matrix <- rawIndex(df, length, cloneCall, coef_matrix)
    }
    coef_matrix$names <- rownames(coef_matrix)
    if (exportTable == TRUE) { return(coef_matrix) }
    coef_matrix <- suppressMessages(melt(coef_matrix))[,-1]
    col <- colorblind_vector(7)
    coef_matrix$variable <- factor(coef_matrix$variable, levels = values)
    coef_matrix$names <- factor(coef_matrix$names, levels = values)
    plot <- ggplot(coef_matrix, aes(x=names, y=variable, fill=value)) +
            geom_tile() + labs(fill = method) +
            geom_text(aes(label = round(value, digits = 3))) +
            scale_fill_gradientn(colors = rev(colorblind_vector(5)), na.value = "white") +
            theme_classic() + theme(axis.title = element_blank())
    return(plot) }
