#' Visualize the number of single cells with clonotype frequencies by cluster
#'
#' View the count of clonotypes frequency group in seurat or SCE object 
#' meta data after combineExpression(). The visualization will take the 
#' new meta data variable "cloneType" and plot the number of cells with
#' each designation using a secondary variable, like cluster. Credit to 
#' the idea goes to Drs. Carmona and Andreatta and their work with
#' \href{https://github.com/carmonalab/ProjecTILs}{ProjectTIL}.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(scRep_example))
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' #Using occupiedscRepertoire()
#' occupiedscRepertoire(sce, x.axis = "ident")
#' table <- occupiedscRepertoire(sce, x.axis = "ident", exportTable = TRUE)
#' 
#' @param sc The single-cell object after \code{\link{combineExpression}}.
#' @param x.axis The variable in the meta data to graph along the x.axis
#' @param label Include the number of clonotype in each category by x.axis variable
#' @param facet.by The column header used for faceting the graph
#' @param proportion Convert the stacked bars into relative proportion
#' @param na.include Visualize NA values or not.
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @importFrom dplyr %>% group_by mutate
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return Stacked bar plot of counts of cells by clonotype frequency group

occupiedscRepertoire <- function(sc, 
                                 x.axis = "ident", 
                                 label = TRUE, 
                                 facet.by = NULL,
                                 proportion = FALSE, 
                                 na.include = FALSE,
                                 exportTable = FALSE, 
                                 palette = "inferno") {
  .checkSingleObject(sc)
  meta <- .grabMeta(sc)
  meta <- melt(table(meta[!is.na(meta$Frequency), 
                          c(x.axis, facet.by, "cloneType")], useNA = "ifany"))
  if (!na.include) {
    meta <- na.omit(meta)
  }
  meta <- meta[meta$value != 0,]
  
  if(proportion) {
    meta <- meta %>%
      group_by(meta[,1]) %>%
      mutate(total = sum(value), 
             prop = value/total)
    meta <- as.data.frame(meta)
  }
  if (exportTable) {
    return(meta)
  }
  col <- length(unique(meta$cloneType))
  if(proportion) {
    plot <- ggplot(meta, aes(x = meta[,x.axis], y = prop, fill = cloneType)) + 
      geom_bar(stat = "identity") 
    lab <- "Proportion of Cells"
    
  } else {
    plot <- ggplot(meta, aes(x = meta[,x.axis], y = value, fill = cloneType)) + 
      geom_bar(stat = "identity") 
    lab <- "Single Cells"
    
  } 
  plot <- plot + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_fill_manual(values = c(.colorizer(palette,col))) + 
    ylab(lab) + 
    theme_classic() + 
    theme(axis.title.x = element_blank())
  if (!is.null(facet.by)) {
    plot <- plot + facet_grid(.~meta[,facet.by])
  }
  if (label) {
    plot <- plot + geom_text(aes(label = value), position = position_stack(vjust = 0.5))
  }
  plot
}