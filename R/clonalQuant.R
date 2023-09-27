#' Quantify the unique clonotypes
#'
#' This function quantifies unique clonotypes. The unique clonotypes 
#' can be either reported as a raw output or scaled to the total number of 
#' clonotypes recovered using the scale parameter. 
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalQuant(combined, cloneCall="strict", scale = TRUE)
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param order Maintain the order of the list when plotting
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2
#' @export
#' @return ggplot of the total or relative unique clonotypes
clonalQuant <- function(df, 
                        cloneCall = "strict", 
                        chain = "both", 
                        scale=FALSE, 
                        group.by = NULL, 
                        split.by = NULL,
                        order = TRUE,
                        exportTable = FALSE, 
                        palette = "inferno") {
  
  if (length(group.by) > 1) { 
    stop("Only one item in the group.by variable can be listed.")
  }
 
  cloneCall <- theCall(cloneCall)
  df <- .data.wrangle(df, split.by, cloneCall, chain)
  
  #Set up mat to store and selecting graph parameters
  if (!is.null(group.by)) {
    x <- group.by
    labs <- group.by
    mat <- data.frame(matrix(NA, length(df), 4))
    colnames(mat) <- c("contigs","values", "total", group.by)
    for (i in seq_along(df)) {
      mat[i,1] <- length(unique(df[[i]][,cloneCall]))
      mat[i,2] <- names(df)[i]
      mat[i,3] <- length(df[[i]][,cloneCall])
      location <- which(colnames(df[[i]]) == group.by)
      mat[i,4] <- df[[i]][1,location] }
    col <- length(unique(mat[,group.by]))
    if (scale) { y <- "scaled"
    mat$scaled <- mat$contigs/mat$total*100
    ylab <- "Percent of Unique Clonotype"
    
    } else { 
      y <- "contigs"
      x <- group.by
      ylab <- "Unique Clonotypes"}
  } else {
    x <- "values"
    labs <- "Samples"
    mat <- data.frame(matrix(NA, length(df), 3))
    colnames(mat) <- c("contigs","values", "total")
    
    #Generating quantification
    for (i in seq_along(df)) {
      mat[i,1] <- length(unique(df[[i]][,cloneCall]))
      mat[i,2] <- names(df)[i]
      mat[i,3] <- length(df[[i]][,cloneCall]) 
    }
    col <- length(unique(mat$values))
    if (scale == TRUE) { 
      y <- "scaled"
      mat$scaled <- mat$contigs/mat$total*100
      ylab <- "Percent of Unique Clonotype"
    } else { 
      y <- "contigs"
      ylab <- "Unique Clonotypes" 
    } 
  }
  if (exportTable) {
    if (length(df) > 1) {
      return(mat)
    }
    # if a single sample, remove the "values" column if NA
    if (is.na(mat[[2]])) {
      mat[[2]] <- NULL
    }
    return(mat)
  }
  
  if(order & is.null(group.by)) {
    mat[,x] <- factor(mat[,x], levels = mat[,x])
  }
  
  plot <- ggplot(aes(x=mat[,x], y=mat[,y], fill=as.factor(mat[,x])), data = mat) +
                  stat_summary(geom = "errorbar", fun.data = mean_se, 
                               position = "dodge", width=.5) + labs(fill = labs) +
                  stat_summary(fun=mean, geom="bar", color="black", lwd=0.25)+
                  theme_classic() + xlab("Samples") + ylab(ylab) +
                  scale_fill_manual(values = .colorizer(palette, col))
  
  # if it is a single run, remove x axis labels if sample name missing
  if ((length(df) == 1) && identical(names(df), NA_character_)) {
    plot <- plot +
      ggplot2::theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  return(plot)
}