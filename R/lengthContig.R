#' Demonstrate the distribution of lengths filtered contigs.
#'
#' This function takes the output of combineTCR(), combineBCR(), or
#' expression2List() and displays either the nucleotide (nt) or amino 
#' acid (aa) sequence length. The sequence length visualized can be 
#' selected using the chains parameter, either the combined clonotype 
#' (both chains) or across all single chains. Visualization can either 
#' be a histogram or if scale = TRUE, the output will be a density plot. 
#' Multiple sequencing runs can be group together using the 
#' group parameter. If a matrix output for the data is preferred, set 
#' exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' lengthContig(combined, cloneCall="aa", chain = "both")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - CDR3 nucleotide (nt), 
#' CDR3 amino acid (aa).
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param order Maintain the order of the list when plotting
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the discrete or relative length distributions of 
#' clonotype sequences
lengthContig <- function(df, 
                         cloneCall = "aa", 
                         chain = "both", 
                         group.by = NULL, 
                         split.by = NULL, 
                         order = TRUE,
                         scale = FALSE, 
                         exportTable = FALSE, 
                         palette = "inferno") {
  df <- list.input.return(df, split.by)
  if(cloneCall == "nt") { cloneCall <- "CTnt"
  ylab <- "CDR3 (NT)"
  } else if (cloneCall == "aa") { cloneCall <- "CTaa" 
  ylab <- "CDR3 (AA)"
  } else { stop("Please make a selection of the type of
                CDR3 sequence to analyze by using `cloneCall`") }
  AB <- grep("TRA", df[[1]]$TCR1)
  IG <- grep("IGH", df[[1]]$TCR1)
  if (length(AB) > 1 && length(IG) == 0) {
    c1 <- "TRA"
    c2 <- "TRB"
  } else if (length(AB) > 1) {
    c1 <- "IGH"
    c2 <- "IGL"
  } else {
    c1 <- "TRG"
    c2 <- "TRD"
  }
  xlab <- "Length"
  Con.df <- NULL
  Con.df <- lengthDF(df, cloneCall, chain, group.by, c1, c2)
  if(is.null(group.by) & order) {
    Con.df[,"values"] <- factor(Con.df[,"values"], levels = names(df))
  }
  names <- names(df)
  if (!is.null(group.by)) { 
    fill <- group.by
    col <- length(unique(Con.df[,group.by]))
    if (scale == TRUE) { yplus <- "Percent of "
    plot <- ggplot(Con.df, aes(fill=Con.df[,group.by],
                               length,(after_stat(count))/sum(after_stat(count))*100)) + 
      geom_density(aes(y=after_stat(scaled)),alpha=.5,lwd=.25,color="black")
    } else { yplus <- "Number of "
    plot<-ggplot(Con.df,aes(as.factor(length),fill=Con.df[,group.by]))+
      geom_bar(position = position_dodge2(preserve = "single"), 
               color="black", lwd=0.25, width=0.9)  +
      scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                                          max(Con.df$length), by = 5),10)) }
  } else if (is.null(group.by)){ 
    fill <- "Samples"
    col <- length(unique(Con.df$values))
    if (scale == TRUE) { yplus <- "Percent of "
    plot <- ggplot(Con.df, aes(length, (after_stat(count))/sum(after_stat(count))*100, 
                               fill=values)) + geom_density(aes(y=after_stat(scaled)), alpha=0.5, 
                                                            lwd=0.25, color="black")
    }  else { yplus <- "Number of "
    plot <- ggplot(Con.df, aes(as.factor(length), fill=values)) +
      geom_bar(position = position_dodge2(preserve = "single"), 
               color="black", lwd=0.25) +
      scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                                          max(Con.df$length), by = 5),10))} }
  plot <- plot + scale_fill_manual(values = .colorizer(palette,col)) +
    labs(fill = fill) + ylab(paste(yplus, ylab, sep="")) +
    xlab(xlab) + theme_classic()
  if (exportTable == TRUE) { return(Con.df) }
  return(plot)}