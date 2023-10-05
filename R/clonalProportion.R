#' Examining the clonal space occupied by specific clonotypes
#'
#' This function calculates the relative clonal space occupied by the 
#' clonotypes. The grouping of these clonotypes is based on the parameter 
#' **clonal.split**, at default, **clonal.split** will group the clonotypes 
#' into bins of 1:10, 11:100, 101:1001, etc. To adjust the clonotypes 
#' selected, change the numbers in the variable split. If a matrix output 
#' for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param clonal.split The cut points for the specific clonotypes.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#'
#' @export
#' @return ggplot of the space occupied by the specific rank of clonotypes
clonalProportion <- function(df,
                             clonal.split = c(10, 100, 1000, 10000, 30000, 100000), 
                             cloneCall = "strict", 
                             chain = "both", 
                             group.by = NULL,
                             exportTable = FALSE, 
                             palette = "inferno") {
    Con.df <- NULL
    cloneCall <- .theCall(cloneCall)
    sco <- is_seurat_object(df) | is_se_object(df)
    df <- .data.wrangle(df, group.by, cloneCall, chain)
    if(!is.null(group.by) & !sco) {
      df <- .groupList(df, group.by)
    }
    
    #Generating data matrix to store value
    mat <- matrix(0, length(df), length(clonal.split), dimnames = list(names(df), 
            paste0('[', c(1, clonal.split[-length(clonal.split)] + 1), ':', clonal.split, ']')))
    
    #Assigning the clonal grouping
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    df <- lapply(df, as.data.frame(table))
    for (i in seq_along(df)) {
        df[[i]] <- rev(sort(as.numeric(df[[i]][,2])))
    }
    cut <- c(1, clonal.split[-length(clonal.split)] + 1)
    for (i in seq_along(clonal.split)) {
        mat[,i] <- vapply(df, function (x) sum(na.omit(x[cut[i]:clonal.split[i]])), 
                            FUN.VALUE = numeric(1))
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    #Plotting
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonal Indices", 
                        values = rev(.colorizer(palette,col))) +
        xlab("Samples") +
        ylab("Occupied Repertoire Space") +
        theme_classic()
    return(plot)
}
