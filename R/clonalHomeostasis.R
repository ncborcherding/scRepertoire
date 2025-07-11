#' Examining the clonal homeostasis of the repertoire
#'
#' This function calculates the space occupied by clone proportions. 
#' The grouping of these clones is based on the parameter **cloneSize**, 
#' at default, **cloneSize** will group the clones into bins of Rare = 0 
#' to 0.0001, Small = 0.0001 to 0.001, etc. To adjust the proportions, 
#' change the number or labeling of the cloneSize parameter. If a matrix 
#' output for the data is preferred, set **exportTable** = TRUE.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalHomeostasis(combined, cloneCall = "gene")
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param cloneSize The cut points of the proportions.
#' @param cloneCall How to call the clone - VDJC gene (**gene**), 
#' CDR3 nucleotide (**nt**), CDR3 amino acid (**aa**),
#' VDJC gene + CDR3 nucleotide (**strict**) or a custom variable 
#' in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#'
#' @export
#' @importFrom utils stack
#' @importFrom stats reshape
#' @concept Visualizing_Clones
#' @return ggplot of the space occupied by the specific proportion of clones
clonalHomeostasis <- function(input.data, 
                              cloneSize = c(Rare = .0001, Small = .001, Medium = .01, Large = .1, Hyperexpanded = 1),
                              cloneCall = "strict", 
                              chain = "both", 
                              group.by = NULL,
                              order.by = NULL,
                              exportTable = FALSE, 
                              palette = "inferno") {
    cloneSize <- c(None = 0, cloneSize)
    input.data <- .dataWrangle(input.data, 
                               group.by, 
                               .theCall(input.data, cloneCall, check.df = FALSE), 
                               chain)
    cloneCall <- .theCall(input.data, cloneCall)
    sco <- .is.seurat.or.se.object(input.data)
    if(!is.null(group.by) & !sco) {
      input.data <- .groupList(input.data, group.by)
    }
    
    #Generating data matrix to store value
    mat <- matrix(0, length(input.data), length(cloneSize) - 1, 
                dimnames = list(names(input.data), 
                names(cloneSize)[-1]))

    #Assigning the clonal grouping
    input.data <- lapply(input.data, '[[', cloneCall)
    input.data <- lapply(input.data, na.omit)
    fun <- function(x) { table(x)/length(x) }
    input.data <- lapply(input.data, fun)
    for (i in 2:length(cloneSize)) {
        mat[,i-1] <- vapply(input.data, function (x) 
                              sum(x[x > cloneSize[i-1] & x <= cloneSize[i]]), FUN.VALUE = numeric(1))
        colnames(mat)[i-1] <- paste0(names(cloneSize[i]), ' (', 
                                    cloneSize[i-1], ' < X <= ', 
                                    cloneSize[i], ')') }
    if (exportTable) { 
      return(mat) 
    }
    
    #Plotting
    mat_df <- as.data.frame(mat) 
    mat_df$Var1 <- rownames(mat_df)
    
    varying_cols <- names(mat_df)[grep("<", names(mat_df))]
    
    mat_melt <- reshape(mat_df,
                        varying = varying_cols,
                        v.names = "value",
                        timevar = "category",
                        times = varying_cols,
                        idvar = "Var1",
                        direction = "long")
    
    if(!is.null(order.by)) {
      mat_melt <- .orderingFunction(vector = order.by,
                                    group.by = "Var1", 
                                    data.frame = mat_melt)
    }
    
    col <- length(unique(mat_melt$category))
    plot <- ggplot(mat_melt, aes(x=as.factor(Var1), y=value, fill=.data[["category"]])) +
        geom_bar(stat = "identity", position="fill", 
                    color = "black", lwd= 0.25) +
        scale_fill_manual(values = rev(.colorizer(palette,col))) +
        labs(x = "Samples", 
             y = "Relative Abundance", 
             fill = "Clonal Group") +
        theme_classic()
    return(plot)
}
