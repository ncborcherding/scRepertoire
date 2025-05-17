#' Quantify the unique clones by group or sample
#'
#' This function quantifies unique clones. The unique clones 
#' can be either reported as a raw output or scaled to the total number of 
#' clones recovered using the scale parameter. 
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalQuant(combined, 
#'             cloneCall="strict", 
#'             scale = TRUE)
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param cloneCall How to call the clone - VDJC gene (**gene**), 
#' CDR3 nucleotide (**nt**), CDR3 amino acid (**aa**),
#' VDJC gene + CDR3 nucleotide (**strict**) or a custom variable 
#' in the data
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param scale Converts the graphs into percentage of unique clones
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @export
#' @concept Visualizing_Clones
#' @return ggplot of the total or relative unique clones
clonalQuant <- function(input.data, 
                        cloneCall = "strict", 
                        chain = "both", 
                        scale=FALSE, 
                        group.by = NULL,
                        order.by = NULL,
                        exportTable = FALSE, 
                        palette = "inferno") {
  
  if (length(group.by) > 1) { 
    stop("Only one item in the group.by variable can be listed.")
  }
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
  
  sco <- is_seurat_object(input.data) | is_se_object(input.data)
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
    mat <- .ordering.function(vector = order.by,
                              group.by = x, 
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
              theme_classic() + xlab("Samples") + 
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
