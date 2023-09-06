#' Quantify the unique clonotypes in the filtered contigs.
#'
#' This function takes the output from combineTCR(), combineBCR(), or 
#' expression2List() and quantifies unique clonotypes. The unique clonotypes 
#' can be either reported as a raw output or scaled to the total number of 
#' clonotypes recovered using the scale parameter. Multiple sequencing 
#' runs can be group together using the group parameter. If a matrix output 
#' for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' quantContig(combined, cloneCall="strict", scale = TRUE)
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param order Maintain the order of the list when plotting
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @import ggplot2
#' @export
#' @return ggplot of the total or relative unique clonotypes
quantContig <- function(df, 
                        cloneCall = "strict", 
                        chain = "both", 
                        scale=FALSE, 
                        group.by = NULL, 
                        split.by = NULL,
                        order = TRUE,
                        exportTable = FALSE, 
                        palette = "inferno") {
  
  if (length(group.by) > 1) { 
    stop("Only one item in the group.by variable can be listed.")
  }
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  if (chain != "both") {
    for(i in seq_along(df)) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  if (!is.null(group.by)) {
    x <- group.by
    labs <- group.by
    Con.df <- data.frame(matrix(NA, length(df), 4))
    colnames(Con.df) <- c("contigs","values", "total", group.by)
    for (i in seq_along(df)) {
      Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
      Con.df[i,2] <- names(df)[i]
      Con.df[i,3] <- length(df[[i]][,cloneCall])
      location <- which(colnames(df[[i]]) == group.by)
      Con.df[i,4] <- df[[i]][1,location] }
    col <- length(unique(Con.df[,group.by]))
    if (scale) { y <- "scaled"
    Con.df$scaled <- Con.df$contigs/Con.df$total*100
    ylab <- "Percent of Unique Clonotype"
    
    } else { y <- "contigs"
    x <- group.by
    ylab <- "Unique Clonotypes"}
  } else {
    x <- "values"
    labs <- "Samples"
    Con.df <- data.frame(matrix(NA, length(df), 3))
    colnames(Con.df) <- c("contigs","values", "total")
    for (i in seq_along(df)) {
      Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
      Con.df[i,2] <- names(df)[i]
      Con.df[i,3] <- length(df[[i]][,cloneCall]) }
    col <- length(unique(Con.df$values))
    if (scale == TRUE) { y <- "scaled"
    Con.df$scaled <- Con.df$contigs/Con.df$total*100
    ylab <- "Percent of Unique Clonotype"
    } else { y <- "contigs"
    ylab <- "Unique Clonotypes" } }
  
  if (exportTable) {
    if (length(df) > 1) {
      return(Con.df)
    }
    
    # if a single sample, remove the "values" column if NA
    if (is.na(Con.df[[2]])) {
      Con.df[[2]] <- NULL
    }
    return(Con.df)
  }
  
  if(order & is.null(group.by)) {
    Con.df[,x] <- factor(Con.df[,x], levels = Con.df[,x])
  }
  plot <- ggplot(aes(x=Con.df[,x], y=Con.df[,y],
                     fill=as.factor(Con.df[,x])), data = Con.df) +
    stat_summary(geom = "errorbar", fun.data = mean_se, 
                 position = "dodge", width=.5) + labs(fill = labs) +
    stat_summary(fun=mean, geom="bar", color="black", lwd=0.25)+
    theme_classic() + xlab("Samples") + ylab(ylab) +
    scale_fill_manual(values = .colorizer(palette, col))
  
  # if it is a single run, remove x axis labels if sample name missing
  if ((length(df) == 1) && identical(names(df), NA_character_)) {
    plot <- plot +
      ggplot2::theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }
  return(plot)
}