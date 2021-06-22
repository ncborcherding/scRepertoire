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
#' @param group The column header for which you would like to analyze the data.
#' @param exportTable Exports a table of the data into the global environment 
#' in addition to the visualization
#' @param n.boots number of bootstraps to downsample in order to get mean diversity
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
#' @author Andrew Malone, Nick Borcherding
clonalDiversity <- function(df, cloneCall = "gene+nt", group = "samples", 
                                exportTable = FALSE, n.boots = 100) {
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  mat <- NULL
  mat_a <- NULL
  sample <- c()
  min <- c()
  for (x in seq_along(df)) {
    min.tmp <- length(which(!is.na(unique(df[[x]][,cloneCall]))))
    min <- c(min.tmp, min)
  }
  min <- min(min)
  if (group == "samples") {
    for (i in seq_along(df)) {
      
      data <- as.data.frame(table(df[[i]][,cloneCall]))
      mat_a <- NULL
      sample <- c()
      
      for (j in seq(seq_len(n.boots))) {
        x<-sample_n(data, min)
        sample <- diversityCall(x)
        mat_a <- rbind.data.frame(mat_a, sample)
      }
      mat_a[is.na(mat_a)] <- 0
      mat_a<- as.data.frame(colMeans(mat_a))
      mat_a<-t(mat_a)
      colnames(mat_a)<-c("a","b","c","d")
      mat <- rbind.data.frame(mat, mat_a)
    }
    
    colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE")
    rownames(mat) <- names(df)
    mat[,group] <- rownames(mat)
    melt <- melt(mat, id.vars = group)
    plot <- ggplot(melt, aes(x = "", y=value)) +
      geom_jitter(shape=21, size=3, width=0.2, aes(fill=melt[,group]))
    
  } else {
    for (i in seq_along(df)) {
      
      data <- as.data.frame(table(df[[i]][,cloneCall]))
      color <- df[[i]][1,group]
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
      mat_a$color <- color
      mat <- rbind(mat, mat_a)
    }
    colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE", group)
    rownames(mat) <- names(df)
    melt <- suppressWarnings(melt(mat, id.vars = group))
    values <- str_sort(as.character(unique(melt[,group])), 
                       numeric = TRUE)
    values2 <- quiet(dput(values))
    melt[,group] <- factor(melt[,group], levels = values2)
    plot <- ggplot(melt, aes(x=melt[,group], y=as.numeric(value))) +
      geom_jitter(shape=21, size=3, width=0.2, 
                  aes(fill=melt[,group])) 
  }
  col <- length(unique(melt[,group]))
  plot <- plot + ylab("Index Score") + scale_fill_manual(name = group, 
                                                         values = colorblind_vector(col)) +
    facet_wrap(~variable, scales = "free", ncol = 4) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (exportTable == TRUE) { return(mat) }
  return(plot) }
