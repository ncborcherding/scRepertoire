#Making lodes to function in alluvial plots
#' @importFrom ggalluvial to_lodes_form
.makingLodes <- function(meta2, color, alpha, facet, set.axes) {
  diffuse_elements <- c()
  if (!is.null(color)) {
    diffuse_elements <- c(diffuse_elements, color)
  }
  if (!is.null(alpha)) {
    diffuse_elements <- c(diffuse_elements, alpha)
  }
  if (!is.null(facet)) {
    diffuse_elements <- c(diffuse_elements, facet)
  }
  
  lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
                         id = "alluvium", axes = set.axes, 
                         diffuse = if(length(diffuse_elements) > 0) diffuse_elements else NULL)
  return(lodes)
}
#' Exploring interaction of clonotypes by single-cell object dynamics
#'
#' View the proportional contribution of clonotypes by seurat or SCE object 
#' meta data after \code{\link{combineExpression}}. The visualization is based on the 
#' ggalluvial package, which requires the aesthetics to be part of the axes 
#' that are visualized. Therefore, alpha, facet, and color should be part of 
#' the the axes you wish to view or will add an additional stratum/column to 
#' the end of the graph.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' #Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident, 1,3)
#' 
#' #Using alluvialClones()
#' alluvialClones(scRep_example, 
#'                    cloneCall = "gene", 
#'                    y.axes = c("Patient", "ident"), 
#'                    color = "ident")
#' 
#' @param sc The single-cell object to visualize after \code{\link{combineExpression}}. 
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt) or CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param y.axes The columns that will separate the proportional 
#' visualizations.
#' @param color The column header or clonotype(s) to be highlighted.
#' @param facet The column label to separate.
#' @param alpha The column header to have gradated opacity.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#' 
#' @import ggplot2
#' @importFrom ggalluvial StatStratum geom_flow geom_stratum to_lodes_form geom_alluvium
#' @importFrom dplyr %>% mutate
#' 
#' @export
#' @return Alluvial ggplot comparing clonotype distribution.
alluvialClones <- function(sc, 
                           cloneCall = "strict", 
                           chain = "both",
                           y.axes = NULL, 
                           color = NULL, 
                           alpha = NULL, 
                           facet = NULL, 
                           exportTable = FALSE,
                           palette = "inferno") {
  
  x <- alluvium <- stratum <- NULL
  .checkSingleObject(sc)
  cloneCall <- .theCall(cloneCall)
  if (length(y.axes) == 0) {
    stop("Make sure you have selected the variable(s) to visualize") 
  }
  meta <- .grabMeta(sc)
  if (chain != "both") {
    meta <- .off.the.chain(meta, chain, cloneCall)
  }
  meta$barcodes <- rownames(meta)
  meta <- meta[!is.na(meta[,cloneCall]),]
  check <- colnames(meta) == color
  if (length(unique(check)) == 1 & unique(check)[1] == FALSE & 
      !is.null(color)) {
    meta <- meta %>% mutate("clonotype(s)" = ifelse(meta[,cloneCall] %in% 
                                                      color, "Selected", "Other"))
    color <- "clonotype(s)" 
  }
  
  #Prepping the data for calculating lodes
  y.axes <- unique(c(y.axes, color, alpha, facet))
  set.axes <- seq_along(y.axes)
  meta2 <- meta[,c(y.axes, cloneCall, "barcodes")]
  meta2 <- unique(na.omit(meta2[!duplicated(as.list(meta2))]))
  
  lodes <- .makingLodes(meta2, color, alpha, facet, set.axes) 
  #Filtering the lodes
  if(any(lodes[,cloneCall] != "")) {
    lodes <- lodes[lodes[,cloneCall] != "",]
  }
  if(any(is.na(lodes[,cloneCall]))) {
    lodes <- lodes[!is.na(lodes[,cloneCall]),]
  }
  if (exportTable) { 
    return(lodes) 
  }
  #Plotting
  plot <- ggplot(data = lodes, aes(x = x, 
                                   stratum = stratum, 
                                   alluvium = alluvium, 
                                   label = stratum)) +
                geom_stratum(width = 0.2) 
    
  if (is.null(color) & is.null(alpha)) {
    plot <- plot + geom_alluvium(width=0.2)
  } else if (!is.null(color) & is.null(alpha)) {
    plot <- plot + geom_flow(aes(fill = lodes[,color]), 
                             stat = "alluvium", 
                             lode.guidance = "forward") + 
                             labs(fill = color,
                             width = 0.2)
  } else if (is.null(color) & !is.null(alpha)) {
    plot <- plot + geom_flow(aes(alpha = lodes[,alpha]), 
                             stat = "alluvium",
                             lode.guidance = "forward") + 
                             labs(alpha = alpha,
                             width = 0.2)
  } else {
    plot <- plot+geom_flow(aes(alpha=lodes[,alpha], 
                               fill=lodes[,color]),
                               stat = "alluvium", 
                               lode.guidance = "forward", 
                               width = 0.2) + 
      labs(fill = color, alpha = alpha) }
  if (length(facet) == 1 & length(facet) < 2) {
    plot <- plot + facet_wrap(.~lodes[,facet], scales="free_y")
  }
  
  plot <- plot + 
            geom_text(stat = ggalluvial::StatStratum, infer.label = FALSE, reverse = TRUE, size = 2) + 
            scale_fill_manual(values = .colorizer(palette,  length(levels(lodes[,color])))) + 
            scale_x_discrete(expand = c(0.025,0.025)) + 
            theme_classic() +
            theme(axis.title.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  line = element_blank())  
  return(plot)
}
