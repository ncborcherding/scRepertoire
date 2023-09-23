#' Demonstrate the relative abundance of clonotypes by group or sample.
#'
#' Displays the number of clonotypes at specific requencies by sample 
#' or group. Visualization can either be a line graph using
#' calculated numbers or if scale = TRUE, the output will be a 
#' density plot. Multiple sequencing runs can be group together 
#' using the group parameter. If a matrix output for the data is 
#' preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' abundanceContig(combined, cloneCall = "gene", scale = FALSE)
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header for which you would like to analyze the data.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param order Maintain the order of the list when plotting
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param exportTable Returns the data frame used for forming the graph
#' to the visualization.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the total or relative abundance of clonotypes 
#' across quanta
abundanceContig <- function(df, 
                            cloneCall = "strict", 
                            chain = "both", 
                            scale=FALSE, 
                            group.by = NULL, 
                            split.by = NULL, 
                            order = TRUE,
                            exportTable = FALSE, 
                            palette = "inferno") {
  df <- list.input.return(df,split.by)
  Con.df <- NULL
  xlab <- "Abundance"
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  names <- names(df)
  for (i in seq_along(df)) {
    if (chain != "both") {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  
  if (!is.null(group.by)) {
    for (i in seq_along(df)) {
      data1 <- parseContigs(df, i, names, cloneCall)
      label <- df[[i]][1,group.by]
      data1[,paste(group.by)] <- label
      Con.df<- rbind.data.frame(Con.df, data1) }
      Con.df <- data.frame(Con.df)
      col <- length(unique(Con.df[,group.by]))
      fill <- group.by
      if (scale == TRUE) { 
        ylab <- "Density of Clonotypes"
        plot <- ggplot(Con.df, aes(x=Abundance, fill=Con.df[,group.by])) +
                      geom_density(aes(y=after_stat(scaled)), alpha=0.5, 
                                   lwd=0.25, color="black", bw=0.5)  +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
        ylab <- "Number of Clonotypes"
        plot <- ggplot(Con.df, aes(x=Abundance, group.by = values, 
                               color = Con.df[,group.by])) +
                        geom_line(stat="count") +
                        scale_color_manual(values = .colorizer(palette,col)) +
                        labs(color = fill)
    }
  } else {
    for (i in seq_along(df)) {
      data1 <- parseContigs(df, i, names, cloneCall)
      Con.df<- rbind.data.frame(Con.df, data1) 
      }
    Con.df <- data.frame(Con.df)
    if(order) {
      Con.df[,"values"] <- factor(Con.df[,"values"], levels = names(df))
    }
    col <- length(unique(Con.df$values))
    fill <- "Samples"
    if (scale == TRUE) { 
      ylab <- "Density of Clonotypes"
      plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                      geom_density(aes(y=after_stat(scaled)), alpha=0.5, lwd=0.25, 
                                   color="black", bw=0.5) +
                      scale_fill_manual(values = .colorizer(palette,col)) +
                      labs(fill = fill)
    } else { 
      ylab <- "Number of Clonotypes"
      plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                               color = values)) +
                      geom_line(stat="count") +
                      scale_color_manual(values = .colorizer(palette,col)) +
                      labs(color = fill)
  } }
  if (exportTable == TRUE) { 
    return(Con.df) 
  }
  plot <- plot + scale_x_log10() + ylab(ylab) + xlab(xlab) +
    theme_classic()
  return(plot) 
}
