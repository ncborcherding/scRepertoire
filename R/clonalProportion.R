#' Examining the clonal space occupied by specific clones
#'
#' This function calculates the relative clonal space occupied by the 
#' clones. The grouping of these clones is based on the parameter 
#' \strong{clonalSplit}, at default, \strong{clonalSplit} will group the clones 
#' into bins of 1:10, 11:100, 101:1001, etc. To adjust the clones 
#' selected, change the numbers in the variable split. If a matrix output 
#' for the data is preferred, set \strong{exportTable} = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param input.data The product of \code{\link{combineTCR}}, 
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}.
#' @param clonalSplit The cut points for the specific clones
#' @param cloneCall How to call the clone - VDJC gene (\strong{gene}), 
#' CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
#' VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable 
#' in the data
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param exportTable Exports a table of the data into the global.
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any 
#' \link[grDevices]{hcl.pals}
#'
#' @import ggplot2
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @importFrom dplyr bind_rows n
#'
#' @export
#' @concept Visualizing_Clones
#' @return ggplot of the space occupied by the specific rank of clones
clonalProportion <- function(input.data,
                             clonalSplit = c(10, 100, 1000, 10000, 30000, 100000), 
                             cloneCall = "strict", 
                             chain = "both", 
                             group.by = NULL,
                             order.by = NULL,
                             exportTable = FALSE, 
                             palette = "inferno") {
    Con.df <- NULL
    input.data <- .data.wrangle(input.data, 
                                group.by, 
                                .theCall(input.data, cloneCall, check.df = FALSE), 
                                chain)
    cloneCall <- .theCall(input.data, cloneCall)
    sco <- is_seurat_object(input.data) | is_se_object(input.data)
    
    if(!is.null(group.by) & !sco) {
      input.data <- .groupList(input.data, group.by)
    }
    
    #Generating data matrix to store value
    mat <- matrix(0, length(input.data), length(clonalSplit), dimnames = list(names(input.data), 
            paste0('[', c(1, clonalSplit[-length(clonalSplit)] + 1), ':', clonalSplit, ']')))
    
    #Assigning the clonal grouping
    input.data <- lapply(input.data, '[[', cloneCall)
    input.data <- lapply(input.data, na.omit)
    input.data <- lapply(input.data, as.data.frame(table))
    for (i in seq_along(input.data)) {
        input.data[[i]] <- rev(sort(as.numeric(input.data[[i]][,2])))
    }
    cut <- c(1, clonalSplit[-length(clonalSplit)] + 1)
    for (i in seq_along(clonalSplit)) {
        mat[,i] <- vapply(input.data, function (x) 
                          sum(na.omit(x[cut[i]:clonalSplit[i]])), FUN.VALUE = numeric(1))
    }
    if (exportTable == TRUE) {
        return(mat)
    }
    #Plotting
    mat_melt <- melt(mat)
    
    if(!is.null(order.by)) {
      mat_melt <- .ordering.function(vector = order.by,
                                     group.by = "Var1", 
                                     data.frame = mat_melt)
    }
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
