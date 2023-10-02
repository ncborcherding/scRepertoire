#' Calculate rarefaction based on the abundance of clones
#' 
#' This functions uses the Hill numbers of order q: species richness (q = 0), 
#' Shannon diversity (q = 1, the exponential of Shannon entropy) and Simpson 
#' diversity (q = 2, the inverse of Simpson concentration) to compute diversity 
#' estimates for rarefaction and extrapolation. The function relies on 
#' \code{\link[iNEXT]{iNEXT}} and read/please cite the 
#' \href{https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12613}{manuscript} 
#' if using this function.
#' 
#' #' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalRarefaction(combined[c(1,2)], cloneCall = "gene", n.boots = 3)
#'
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The variable to use for grouping
#' @param plot.type sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).   
#' @param hill.numbers The Hill numbers to be plotted out (0 - species richness, 1 - Shannon, 2 - Simpson)
#' @param n.boots The number of bootstraps to downsample in order to get mean diversity
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#' 
#' @importFrom iNEXT iNEXT ggiNEXT
#' @import ggplot2
#' @export
clonalRarefaction <- function(df,
                              cloneCall = "strict", 
                              chain = "both", 
                              group.by = NULL, 
                              plot.type = 1,
                              hill.numbers = 0,
                              n.boots = 20,
                              exportTable = FALSE,
                              palette = "inferno") {
  cloneCall <- .theCall(cloneCall)
  
  sco <- is_seurat_object(df) | is_se_object(df)
  df <- .data.wrangle(df, group.by, "CTgene", chain)
  if(!is.null(group.by) & !sco) {
    df <- .groupList(df, group.by)
  }
  
  mat.list <- lapply(df, function(x) {
                  table(x[,cloneCall])
  })
  col <- length(df)
  mat <- iNEXT(mat.list, q=hill.numbers, datatype="abundance",nboot = n.boots) 
  plot <- ggiNEXT(mat, type=plot.type)  + 
            scale_shape_manual(values = rep(16,col)) + 
            scale_fill_manual(values = rep("white", col)) + 
            scale_color_manual(values = c(.colorizer(palette,col)))  + 
            theme_classic()
  if (exportTable == TRUE) { 
    return(plot["data"]) 
  } else {
    return(plot)
  }
  
}
