#' Visualize distribution of clonal frequency overlaid on dimensional reduction plots
#'
#' This function allows the user to visualize the clonal expansion by overlaying the 
#' cells with specific clonal frequency onto the dimensional reduction plots in Seurat.
#' Credit to the idea goes to Drs Andreatta and Carmona and their work with
#' \href{https://github.com/carmonalab/ProjecTILs}{ProjectTIL}.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' #Using combineExpresion()
#' scRep_example <- combineExpression(combined, 
#'                                    scRep_example)
#' 
#' #Using clonalOverlay()
#' clonalOverlay(scRep_example, 
#'               reduction = "umap", 
#'               freq.cutpoint = 0.3, 
#'               bins = 5) 
#' 
#' @param sc The single-cell object after \code{\link{combineExpression}}.
#' @param reduction The dimensional reduction to visualize
#' @param freq.cutpoint The overlay cutpoint to include, this corresponds to the 
#' Frequency variable in the single-cell objecter
#' @param bins The number of contours to the overlay
#' @param facet meta data variable to facet the comparison
#' 
#' @import ggplot2
#' @importFrom SeuratObject Embeddings
#' @export
#' @author Francesco Mazziotta, Nick Borcherding
#' 
#' @return ggplot object

clonalOverlay <- function(sc, 
                          reduction = NULL, 
                          freq.cutpoint = 30, 
                          bins = 25, 
                          facet = NULL) {
  checkSingleObject(sc)
  tmp <- data.frame(grabMeta(sc), identity = sc@active.ident, get.coord(sc, reduction))
  if (!is.null(facet)) {
    facet <- tmp[,facet]
    tmp <- data.frame(facet, tmp)
  }
  tmp$include <- ifelse(tmp$Frequency >= freq.cutpoint, "Yes", NA)
  tmp2 <- subset(tmp, include == "Yes")
  plot <- ggplot(tmp2, mapping = aes(x=tmp2[,(ncol(tmp2)-2)], y = tmp2[,(ncol(tmp2)-1)])) +
    geom_point(tmp, mapping = aes(x=as.numeric(tmp[,(ncol(tmp)-2)]), y = as.numeric(tmp[,(ncol(tmp)-1)]), color = tmp[,"identity"]), size= 0.5) +
    geom_density_2d(color = "black", lwd=0.25, bins = bins) + 
    theme_classic() +
    labs(color = "Active Identity") +
    xlab("Dimension 1") + ylab("Dimension 2")
  if (!is.null(facet)) {
    plot <- plot + facet_wrap(~facet) 
  }
  plot
}