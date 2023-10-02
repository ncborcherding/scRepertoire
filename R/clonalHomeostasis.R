#' Examining the clonal homeostasis
#'
#' This function calculates the space occupied by clonotype proportions. 
#' The grouping of these clonotypes is based on the parameter cloneSize, 
#' at default, cloneSize will group the clonotypes into bins of Rare = 0 
#' to 0.0001, Small = 0.0001 to 0.001, etc. To adjust the proportions, 
#' change the number or labeling of the cloneSize paramter. If a matrix 
#' output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalHomeostasis(combined, cloneCall = "gene")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneSize The cutpoints of the proportions.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows
#' @export
#' @return ggplot of the space occupied by the specific proportion of clonotypes
clonalHomeostasis <- function(df, 
                              cloneSize = 
                                c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1),
                              cloneCall = "strict", 
                              chain = "both", 
                              group.by = NULL,
                              exportTable = FALSE, 
                              palette = "inferno") {
    cloneSize <- c(None = 0, cloneSize)
    
    cloneCall <- .theCall(cloneCall)
    sco <- is_seurat_object(df) | is_se_object(df)
    df <- .data.wrangle(df, group.by, cloneCall, chain)
    if(!is.null(group.by) & !sco) {
      df <- .groupList(df, group.by)
    }
    
    #Generating data matrix to store value
    mat <- matrix(0, length(df), length(cloneSize) - 1, 
                dimnames = list(names(df), 
                names(cloneSize)[-1]))

    #Assigning the clonal grouping
    df <- lapply(df, '[[', cloneCall)
    df <- lapply(df, na.omit)
    fun <- function(x) { table(x)/length(x) }
    df <- lapply(df, fun)
    for (i in 2:length(cloneSize)) {
        mat[,i-1] <- vapply(df, function (x) sum(x[x > cloneSize[i-1] & x <= 
                            cloneSize[i]]), FUN.VALUE = numeric(1))
        colnames(mat)[i-1] <- paste0(names(cloneSize[i]), ' (', 
                                    cloneSize[i-1], ' < X <= ', 
                                    cloneSize[i], ')') }
    if (exportTable) { 
      return(mat) 
    }
    
    #Plotting
    mat_melt <- melt(mat)
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=Var2)) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(name = "Clonotype Group", 
                    values = rev(.colorizer(palette,col))) +
        xlab("Samples") +
        ylab("Relative Abundance") +
        theme_classic()
    return(plot)
}
