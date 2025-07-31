#' Calculate rarefaction based on the abundance of clones
#' 
#' This functions uses the Hill numbers of order q: species richness (`q = 0`), 
#' Shannon diversity (`q = 1`), the exponential of Shannon entropy and Simpson 
#' diversity (`q = 2`, the inverse of Simpson concentration) to compute diversity 
#' estimates for rarefaction and extrapolation. The function relies on the
#' [iNEXT::iNEXT()] R package. Please read and cite the 
#' [manuscript](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12613) 
#' if using this function. The input into the iNEXT calculation is abundance, 
#' incidence-based calculations are not supported.
#' 
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#'                                     
#' # Using clonalRarefaction()
#' clonalRarefaction(combined[c(1,2)], cloneCall = "gene", n.boots = 3)
#'
#'
#' @param input.data The product of [combineTCR()], [combineBCR()], or 
#' [combineExpression()].
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param plot.type sample-size-based rarefaction/extrapolation curve 
#' (`type = 1`); sample completeness curve (`type = 2`); 
#' coverage-based rarefaction/extrapolation curve (`type = 3`).   
#' @param hill.numbers The Hill numbers to be plotted out 
#' (0 - species richness, 1 - Shannon, 2 - Simpson)
#' @param n.boots The number of bootstrap replicates used to derive confidence 
#' intervals for the diversity estimates. More replicates can provide a more 
#' reliable measure of statistical variability.
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#' @param ... Additional arguments passed to the ggplot theme
#'
#' @importFrom iNEXT iNEXT ggiNEXT
#' @export
#' @concept Visualizing_Clones
clonalRarefaction <- function(input.data,
                              cloneCall = "strict", 
                              chain = "both", 
                              group.by = NULL, 
                              plot.type = 1,
                              hill.numbers = 0,
                              n.boots = 20,
                              exportTable = FALSE,
                              palette = "inferno",
                              ...) {
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
  
  mat.list <- lapply(input.data, function(x) {
                  table(x[,cloneCall])
  })
  col <- length(input.data)
  
  # Compute diversity estimates along with bootstrap-derived confidence intervals.
  # The iNEXT function uses 'n.boots' replicates to estimate variability.
  mat <- iNEXT(mat.list, q=hill.numbers, datatype="abundance",nboot = n.boots) 
  plot <- suppressMessages(ggiNEXT(mat, type=plot.type) + 
            scale_shape_manual(values = rep(16,col)) + 
            scale_fill_manual(values = c(.colorizer(palette,col))) + 
            scale_color_manual(values = c(.colorizer(palette,col)))  + 
            .themeRepertoire(...))
  
  if (exportTable) { 
    return(plot["data"]) 
  } else {
    return(plot)
  }
  
}
