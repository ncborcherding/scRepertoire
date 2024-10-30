#' Visualize the number of single cells with cloneSizes by cluster
#'
#' View the count of clones frequency group in Seurat or SCE object 
#' meta data after [combineExpression()]. The visualization 
#' will take the new meta data variable **"cloneSize"** and 
#' plot the number of cells with each designation using a secondary 
#' variable, like cluster. Credit to the idea goes to Drs. Carmona 
#' and Andreatta and their work with [ProjectTIL](https://github.com/carmonalab/ProjecTILs).
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
#' 
#' #Using clonalOccupy
#' clonalOccupy(scRep_example, x.axis = "ident")
#' table <- clonalOccupy(scRep_example, x.axis = "ident", exportTable = TRUE)
#' 
#' @param sc.data The single-cell object after [combineExpression()]
#' @param x.axis The variable in the meta data to graph along the x.axis.
#' @param label Include the number of clone in each category by x.axis variable
#' @param facet.by The column header used for faceting the graph
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order description
#' @param proportion Convert the stacked bars into relative proportion
#' @param na.include Visualize NA values or not
#' @param exportTable Exports a table of the data into the global 
#' environment in addition to the visualization
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @importFrom dplyr %>% group_by mutate count
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @concept SC_Functions
#' @return Stacked bar plot of counts of cells by clone frequency group

clonalOccupy <- function(sc.data, 
                         x.axis = "ident", 
                         label = TRUE, 
                         facet.by = NULL,
                         order.by = NULL,
                         proportion = FALSE, 
                         na.include = FALSE,
                         exportTable = FALSE, 
                         palette = "inferno") {
  .checkSingleObject(sc.data)
  meta <- .grabMeta(sc.data)
  
  meta <- meta %>%
            group_by(meta[,x.axis], meta[,facet.by], cloneSize) %>%
            count() %>%
            as.data.frame()
  colnames(meta)[1] <- x.axis
  
  if(!is.null(order.by)) {
    meta <- .ordering.function(vector = order.by,
                               group.by = x.axis, 
                               data.frame = meta)
  } else {
    meta[,1] <- as.factor(meta[,1])
  }
  
  colnames(meta)[1] <- x.axis
  if(!is.null(facet.by)) {
    colnames(meta)[2] <- facet.by
  }
 
  #Check for NAs
  if (!na.include) {
    meta <- na.omit(meta)
  }
  meta <- meta[meta$n != 0,]
  
  #Convert to proportion
  if(proportion) {
    meta <- meta %>%
      group_by(meta[,1]) %>%
      mutate(total = sum(n), 
             prop = n/total)
    meta <- as.data.frame(meta)
  }
  if (exportTable) {
    return(meta)
  }
  #Plotting
  col <- length(unique(meta[,"cloneSize"]))
  if(proportion) {
    plot <- ggplot(meta, aes(x = meta[,x.axis], y = prop, fill = cloneSize)) + 
              geom_bar(stat = "identity", color = "black", lwd = 0.25) 
    lab <- "Proportion of Cells"
    
  } else {
    plot <- ggplot(meta, aes(x = meta[,x.axis], y = n, fill = cloneSize)) + 
              geom_bar(stat = "identity", color = "black", lwd = 0.25) 
    lab <- "Single Cells"
    
  } 
  plot <- plot + 
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
            scale_fill_manual(values = rev(c(.colorizer(palette,col)))) + 
            ylab(lab) + 
            theme_classic() + 
            theme(axis.title.x = element_blank())
  if (!is.null(facet.by)) {
    plot <- plot + 
              facet_grid(.~meta[,facet.by])
  }
  if (label) {
    plot <- plot + 
              geom_text(aes(label = n), position = position_stack(vjust = 0.5))
  }
  plot
}
