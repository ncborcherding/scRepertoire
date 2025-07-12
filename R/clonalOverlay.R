#' Visualize distribution of clonal frequency overlaid on dimensional reduction plots
#'
#' This function allows the user to visualize the clonal expansion by overlaying the 
#' cells with specific clonal frequency onto the dimensional reduction plots in Seurat.
#' Credit to the idea goes to Drs Andreatta and Carmona and their work with
#' [ProjectTIL](https://github.com/carmonalab/ProjecTILs).
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' # Using combineExpresion()
#' scRep_example <- combineExpression(combined, 
#'                                    scRep_example)
#' 
#' # Using clonalOverlay()
#' clonalOverlay(scRep_example, 
#'               reduction = "umap", 
#'               cutpoint = 3, 
#'               bins = 5) 
#' 
#' @param sc.data The single-cell object after [combineExpression()].
#' @param reduction The dimensional reduction to visualize.
#' @param cut.category Meta data variable of the single-cell object to use for 
#' filtering.
#' @param cutpoint The overlay cut point to include, this corresponds to the 
#' cut.category variable in the meta data of the single-cell object.
#' @param bins The number of contours to the overlay
#' @param pt.size The point size for plotting (default is 0.5)
#' @param pt.alpha The alpha value for plotting (default is 1)
#' @param facet.by meta data variable to facet the comparison
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @importFrom SeuratObject Embeddings
#' @export
#' @concept SC_Functions
#' @author Francesco Mazziotta, Nick Borcherding
#' 
#' @return ggplot object

clonalOverlay <- function(sc.data, 
                          reduction = NULL, 
                          cut.category = "clonalFrequency",
                          cutpoint = 30, 
                          bins = 25, 
                          pt.size = 0.5,
                          pt.alpha = 1,
                          facet.by = NULL,
                          ...) {
  .checkSingleObject(sc.data)

  #Forming the data frame to plot
  tmp <- data.frame(.grabMeta(sc.data), .getCoord(sc.data, reduction))
  
  if(!cut.category %in% colnames(tmp)) {
    stop("If filtering the data using a cutpoint, ensure the cut.category correspond to a variable in the meta data.")
  }
  #Add facet variable if present
  if (!is.null(facet.by)) { 
    facet.by <- tmp[,facet.by]
    tmp <- data.frame(facet.by, tmp)
  }
  #If using cut.category for filtering
  if(!is.null(cut.category) & !is.null(cutpoint)) {
    tmp$include <- ifelse(tmp[,cut.category] >= cutpoint, "Yes", NA)
    tmp2 <- subset(tmp, include == "Yes")
  }
  
  #Plotting
  plot <- ggplot(tmp2, mapping = aes(x = tmp2[,(ncol(tmp2)-2)], 
                                     y = tmp2[,(ncol(tmp2)-1)])) +
    geom_point(tmp, mapping = aes(x = as.numeric(tmp[,(ncol(tmp)-2)]), 
                                  y = as.numeric(tmp[,(ncol(tmp)-1)]), 
                                  color = tmp[,"ident"]), 
               size = pt.size,
               alpha = pt.alpha) +
    geom_density_2d(color = "black", lwd=0.25, bins = bins) + 
    .themeRepertoire() + 
    labs(color = "Active Identity") +
    xlab("Dimension 1") + 
    ylab("Dimension 2")
  if (!is.null(facet.by)) {
    plot <- plot + facet_wrap(~facet.by) 
  }
  return(plot)
}
