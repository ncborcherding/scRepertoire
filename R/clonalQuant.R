#' Plot Number or Proportions of Clones
#'
#' This function quantifies unique clones. The unique clones 
#' can be either reported as a raw output or scaled to the total number of 
#' clones recovered using the scale parameter. 
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using clonalQuant()
#' clonalQuant(combined, 
#'             cloneCall="strict", 
#'             scale = TRUE)
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param scale Converts the graphs into percentage of unique clones
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @param ... Additional arguments passed to the ggplot theme
#' 
#' @export
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing the total or relative number of clones 
#' or a data.frame if `exportTable = TRUE`.
clonalQuant <- function(input.data, 
                        cloneCall = "strict", 
                        chain = "both", 
                        scale=FALSE, 
                        group.by = NULL,
                        order.by = NULL,
                        exportTable = FALSE, 
                        palette = "inferno",
                        ...) {
  
  if (length(group.by) > 1) { 
    stop("Only one item in the group.by variable can be listed.")
  }
  input.data <- .dataWrangle(input.data, 
                             group.by, 
                             .theCall(input.data, cloneCall, check.df = FALSE), 
                             chain)
  cloneCall <- .theCall(input.data, cloneCall)
  
  sco <- .is.seurat.or.se.object(input.data)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  mat.names <- c("contigs","values", "total", group.by)
  #Set up mat to store and selecting graph parameters
  if (!is.null(group.by)) {
    x <- group.by
    labs <- group.by
  } else {
    x <- "values"
    labs <- "Samples"
    col <- length(unique(names(input.data)))
  }
  mat <- data.frame(matrix(NA, length(input.data), length(mat.names)))
  colnames(mat) <- mat.names
  for (i in seq_along(input.data)) {
      mat[i,1] <- length(na.omit(unique(input.data[[i]][,cloneCall])))
      mat[i,2] <- names(input.data)[i]
      mat[i,3] <- length(na.omit(input.data[[i]][,cloneCall]))
      if (!is.null(group.by)) {
        location <- which(colnames(input.data[[i]]) == group.by)
        mat[i,4] <- as.vector(input.data[[i]][1,location])
      }
  }
  if (scale) { 
      y <- "scaled"
      mat$scaled <- mat$contigs/mat$total*100
      ylab <- "Percent of Unique Clones"
   } else { 
      y <- "contigs"
      ylab <- "Unique Clones"
   }
  
  if (exportTable) {
    if (length(input.data) > 1) {
      return(mat)
    }
    # if a single sample, remove the "values" column if NA
    if (is.na(mat[[2]])) {
      mat[[2]] <- NULL
    }
    return(mat)
  }
  
  if(!is.null(group.by)) {
    col <- length(unique(mat[,group.by]))
  }
  mat[,x] = factor(mat[,x], levels = names(input.data))
  
  if(!is.null(order.by)) {
    mat <- .orderingFunction(vector = order.by,
                             group.by = "values", 
                             data.frame = mat)
  }
  
  #Plotting
  plot <- ggplot(data = mat, 
                 aes(x=mat[,x], y=mat[,y], fill=mat[,x])) +
            stat_summary(geom = "errorbar", 
                         fun.data = mean_se, 
                         position = "dodge", 
                         width=.5) + 
              labs(fill = labs) +
              ylab(ylab) +
              stat_summary(fun=mean, geom="bar", color="black", lwd=0.25)+
              .themeRepertoire(...) + 
              xlab("Samples") + 
              scale_fill_manual(values = .colorizer(palette, col))
  
  # if it is a single run, remove x axis labels if sample name missing
  if ((length(input.data) == 1) && identical(names(input.data), NA_character_)) {
    plot <- plot +
      ggplot2::theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  return(plot)
}
