#' Examine the clonal diversity of samples
#'
#' This function calculates traditional measures of diversity - Shannon, 
#' inverse Simpson, Chao1 index, and abundance-based coverage estimators 
#' (ACE) by sample or group. The function automatically down samples the
#' diversity metrics using 100 boot straps The group paramter can be 
#' used to condense the individual samples. If a matrix output for 
#' the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalDiversity(combined, cloneCall = "gene")
#'
#' @param df The product of combineTCR(), combineBCR(),  or expression2List().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param groupBy Variable in which to group the diversity calculation
#' @param x.axis Additional variable in which to split the x.axis
#' @param exportTable Exports a table of the data into the global environment 
#' in addition to the visualization
#' @param n.boots number of bootstraps to downsample in order to get mean diversity
#' @importFrom stringr str_sort str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
#' @author Andrew Malone, Nick Borcherding
clonalDiversity <- function(df, cloneCall = "gene+nt", chain = "both",
                            groupBy = NULL, x.axis = NULL,
                            exportTable = FALSE, n.boots = 100) {
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  if (chain != "both") {
    df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
  }
  min <- c()
  mat <- NULL
  mat_a <- NULL
  sample <- c()
  if (!is.null(groupBy) || !is.null(x.axis)) {
    df <- bind_rows(df, .id = "element.names")
    df$group.element <- paste0(df[,groupBy], ".", df[,x.axis])
    group.element.uniq <- unique(df$group.element)
    df <- split(df, f = df[,"group.element"])
  }
  for (x in seq_along(df)) {
    min.tmp <- length(which(!is.na(unique(df[[x]][,cloneCall]))))
    min <- c(min.tmp, min)
  }
  min <- min(min)
  
  for (i in seq_along(df)) {
      data <- as.data.frame(table(df[[i]][,cloneCall]))
      mat_a <- NULL
      sample <- c()
      
      for (j in seq(seq_len(n.boots))) {
        x <- sample_n(data, min)
        sample <- diversityCall(x)
        mat_a <- rbind(mat_a, sample)
      }
      mat_a[is.na(mat_a)] <- 0
      mat_a<- colMeans(mat_a)
      mat_a<-as.data.frame(t(mat_a))
      mat <- rbind(mat, mat_a)
    }
    colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE")
    if (!is.null(groupBy)) {
      mat[,groupBy] <- str_split(group.element.uniq, "[.]", simplify = TRUE)[,1]
    } else {
      groupBy <- "Group"
      mat[,groupBy] <- names(df)
    }
    if (!is.null(x.axis)) {
      mat[,x.axis] <- str_split(group.element.uniq, "[.]", simplify = TRUE)[,2]
    } else {
      x.axis <- "x.axis"
      mat[,x.axis] <- 1
    }
    rownames(mat) <- names(df)
  
    melt <- suppressMessages(melt(mat, id.vars = c(groupBy, x.axis)))
    values <- str_sort(as.character(unique(melt[,groupBy])), 
                       numeric = TRUE)
    values <- quiet(dput(values))
    melt[,groupBy] <- factor(melt[,groupBy], levels = values)
    if (x.axis == "x.axis") {
        plot <- ggplot(melt, aes(x=1, y=as.numeric(value)))
    } else {
      plot <- ggplot(melt, aes(x=melt[,x.axis], y=as.numeric(value)))
    }
    plot <- plot +
      geom_boxplot() +
      geom_jitter(aes(color = melt[,groupBy]), size = 3) + 
      labs(color="Group") +
      ylab("Index Score") +
      scale_color_manual(values = colorblind_vector(length(unique(melt[,groupBy])))) +
    facet_wrap(~variable, scales = "free", ncol = 4) +
      theme_classic() + 
      theme(axis.title.x = element_blank())
    if (x.axis == "x.axis") { 
      plot <- plot + theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
      }
  if (exportTable == TRUE) { return(mat) }
  return(plot) 
}
