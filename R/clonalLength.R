#' Demonstrate the distribution of clonal length
#'
#' This function displays either the nucleotide (**nt**) or amino 
#' acid (**aa**) sequence length. The sequence length visualized 
#' can be selected using the chains parameter, either the combined clone 
#' (both chains) or across all single chains. Visualization can either 
#' be a histogram or if **scale** = TRUE, the output will 
#' be a density plot. Multiple sequencing runs can be group together 
#' using the group.by parameter.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalLength(combined, cloneCall="aa", chain = "both")
#'
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()]
#' @param cloneCall How to call the clone - CDR3 nucleotide (**nt**) 
#' or CDR3 amino acid (**aa**)
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order description
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals]
#' @export
#' @concept Visualizing_Clones
#' @return ggplot of the discrete or relative length distributions of
#' clone sequences
clonalLength <- function(input.data, 
                         cloneCall = "aa", 
                         chain = "both", 
                         group.by = NULL, 
                         order.by = NULL,
                         scale = FALSE, 
                         exportTable = FALSE, 
                         palette = "inferno") {
  
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
  
  #Sorting out graphing parameters
  xlab <- "Length"
  if(cloneCall == "CTnt") { 
    ylab <- "CDR3 (NT)"
  } else if (cloneCall == "CTaa") { 
    ylab <- "CDR3 (AA)"
  } else { 
    stop("Please make a selection of the type of
          CDR3 sequence to analyze by using `cloneCall`")
  }
  
  #Calculating Length
  Con.df <- NULL
  Con.df <- .lengthDF(input.data, cloneCall, chain, group.by)
  
  names <- names(input.data)
  
  #Skip plotting if want to export table
  if (exportTable == TRUE) { 
    return(Con.df) 
  }
  
  if(!is.null(order.by)) {
    if (!is.null(group.by)) { 
      Con.df <- .ordering.function(vector = order.by,
                                   group.by = group.by, 
                                   data.frame = Con.df)
    } else {
      Con.df <- .ordering.function(vector = order.by,
                                   group.by = "values", 
                                   data.frame = Con.df)
    }
  }
  
  #Plotting
  if (!is.null(group.by)) { 
    fill <- group.by
    col <- length(unique(Con.df[,group.by]))
    if (scale == TRUE) { 
      yplus <- "Percent of "
      plot <- ggplot(Con.df, aes(fill=Con.df[,group.by],
                                 x = length,
                                 y = (after_stat(count))/sum(after_stat(count))*100)) + 
                     geom_density(aes(y=after_stat(scaled)),
                                  alpha=.5, lwd=.25, color="black")
    } else { 
      yplus <- "Number of "
      plot <- ggplot(Con.df,aes(x = as.factor(length),
                                fill=Con.df[,group.by]))+
                     geom_bar(position = position_dodge2(preserve = "single"), 
                              color="black", lwd=0.25, width=0.9)  +
                      scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                                          max(Con.df$length), by = 5),10)) }
  } else if (is.null(group.by)){ 
      fill <- "Samples"
      col <- length(unique(Con.df$values))
      if (scale == TRUE) { 
        yplus <- "Density of "
        plot <- ggplot(Con.df, aes(x = length, 
                                   y = (after_stat(count))/sum(after_stat(count))*100, 
                                   fill=values)) + 
                        geom_density(aes(y=after_stat(scaled)), 
                                     alpha=0.5, lwd=0.25, color="black")
      }  else { 
        yplus <- "Number of "
        plot <- ggplot(Con.df, aes(as.factor(length), fill=values)) +
                        geom_bar(position = position_dodge2(preserve = "single"), 
                                 color="black", lwd=0.25) +
                        scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                                          max(Con.df$length), by = 5),10))} }
  plot <- plot + 
          scale_fill_manual(values = .colorizer(palette,col)) +
          labs(fill = fill) + 
          ylab(paste(yplus, ylab, sep="")) +
          xlab(xlab) + 
          theme_classic()
  
  return(plot)
}
