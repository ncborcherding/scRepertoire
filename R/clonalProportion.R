#' Plot the Clonal Space Occupied by Specific Clones
#'
#' This function calculates the relative clonal space occupied by the 
#' clones. The grouping of these clones is based on the parameter 
#' `clonalSplit`, at default, `clonalSplit` will group the clones 
#' into bins of 1:10, 11:100, 101:1001, etc. To adjust the clones 
#' selected, change the numbers in the variable split. If a matrix output 
#' for the data is preferred, set `exportTable` = TRUE.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' # Using clonalProportion()
#' clonalProportion(combined, cloneCall = "gene")
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param clonalSplit The cut points for the specific clones, default = c(10, 
#' 100, 1000, 10000, 30000, 1e+05)
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @export
#' @concept Visualizing_Clones
#' @return A ggplot object dividing space occupied by ranks of clones or a 
#' data.frame if `exportTable = TRUE`.
clonalProportion <- function(input.data,
                             clonalSplit = c(10, 100, 1000, 10000, 30000, 100000), 
                             cloneCall = "strict", 
                             chain = "both", 
                             group.by = NULL,
                             order.by = NULL,
                             exportTable = FALSE, 
                             palette = "inferno",
                             ...) {
    Con.df <- NULL
    input.data <- .dataWrangle(input.data, 
                               group.by, 
                               .theCall(input.data, cloneCall, 
                                        check.df = FALSE, silent = TRUE), 
                               chain)
    cloneCall <- .theCall(input.data, cloneCall)
    sco <- .is.seurat.or.se.object(input.data)
    
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
    # Plotting
    mat_melt <- data.frame(
      Var1 = rep(rownames(mat), ncol(mat)),
      Var2 = rep(colnames(mat), each = nrow(mat)),
      value = as.vector(mat)
    )
    mat_melt[["Var2"]] <- factor(mat_melt[["Var2"]], 
                                 levels = colnames(mat))
    
    if(!is.null(order.by)) {
      mat_melt <- .orderingFunction(vector = order.by,
                                    group.by = "Var1", 
                                    data.frame = mat_melt)
    }
    col <- length(unique(mat_melt$Var2))
    plot <- ggplot(mat_melt, aes(x=as.factor(.data[["Var1"]]), y=.data[["value"]], fill=.data[["Var2"]])) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(values = rev(.colorizer(palette,col))) +
        labs(x = "Samples",
             y = "Occupied Repertoire Space", 
             fill = "Clonal Indices") +
      .themeRepertoire(...)
    return(plot)
}
