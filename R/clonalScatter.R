#' Scatter plot comparing the expansion of two samples
#'
#' This function produces a scatter plot directly comparing the specific clonotypes
#' between two samples. The clonotypes will be categorized by counts into singlets or multiplets, 
#' either exclusive or shared between the selected samples. Visualization inspired 
#' by the work of \href{https://pubmed.ncbi.nlm.nih.gov/32103181/}{Wu, T, et al 2020}. 
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalScatter(combined, 
#'               x.axis = "P17B", 
#'               y.axis = "P17L",
#'               graph = "proportion")
#' 
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param x.axis name of the list element to appear on the x.axis
#' @param y.axis name of the list element to appear on the y.axis
#' @param dot.size either total or the name of the list element to 
#' use for size of dots
#' @param group.by The variable to use for grouping
#' @param graph graph either proportion or raw clonotype count
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' 
#' @import ggplot2
#' 
#' @export
#' @return ggplot of the relative clonotype numbers

clonalScatter <- function(df, 
                          cloneCall ="strict", 
                          x.axis = NULL, 
                          y.axis = NULL,
                          chain = "both",
                          dot.size = "total", 
                          group.by = NULL,
                          graph = "proportion", 
                          exportTable = FALSE,
                          palette = "inferno") {
  cloneCall <- .theCall(cloneCall)
  df <- .data.wrangle(df, group.by, cloneCall, chain)
  axes <- which(names(df) %in% c(x.axis, y.axis, dot.size))
  
  #Making new data frame and adding variables to graph
  x.df <- as.data.frame(table(df[[x.axis]][,cloneCall]))
  colnames(x.df)[2] <- x.axis
  y.df <- as.data.frame(table(df[[y.axis]][,cloneCall]))
  colnames(y.df)[2] <- y.axis
  mat <- merge(x.df, y.df, by = "Var1", all = TRUE)
  
  if (dot.size != "total") {
    if (dot.size %!in% colnames(mat)) {
      size.df <- as.data.frame(table(df[[dot.size]][,cloneCall]))
      colnames(size.df)[2] <- dot.size
      mat <- merge(mat, size.df, by = "Var1", all = TRUE) }
    mat[is.na(mat)] <- 0
    mat[,paste0("size", ".fraction")] <- mat[,dot.size]/sum(mat[,dot.size])
    labeling <- unique(c(x.axis, y.axis, dot.size))
  } else {
    mat[is.na(mat)] <- 0
    labeling <- unique(c(x.axis, y.axis)) }
  
  mat[,"class"] <- NA
  mat[,"sum"] <- rowSums(mat[,labeling])
  
  #Assigning class based on x and y-axis counts
  for (i in seq_along(labeling)) {
    if (length(labeling) > 2) {
      mat[,"class"] <- ifelse(mat[,labeling[i]] == 1 & rowSums(mat[,labeling[-i]]) == 0, 
                                      paste0(labeling[i], ".singlet"), #if
                                      mat[,"class"]) #else
      mat[,"class"] <- ifelse(mat[,labeling[i]] > 1 & rowSums(mat[,labeling[-i]]) == 0, 
                                      paste0(labeling[i], ".expanded"), #if
                                      mat[,"class"]) #else
    } else if (length(labeling) == 2) {
      mat[,"class"] <- ifelse(mat[,labeling[i]] == 1 & mat[,labeling[-i]] == 0, 
                                      paste0(labeling[i], ".singlet"), #if
                                      mat[,"class"]) #else
      mat[,"class"] <- ifelse(mat[,labeling[i]] > 1 & mat[,labeling[-i]] == 0, 
                                      paste0(labeling[i], ".expanded"), #if
                                      mat[,"class"]) #else
    }
  }
  #Adding dual-expanded class
  mat[,"class"] <- ifelse(mat[,y.axis] >= 1 & mat[,x.axis] >= 1, paste0("dual.expanded"), mat[,"class"])
  
  #Calculating relative proportion
  mat[,paste0(x.axis, ".fraction")] <- mat[,x.axis]/sum(mat[,x.axis])
  mat[,paste0(y.axis, ".fraction")] <- mat[,y.axis]/sum(mat[,y.axis])
  
  #Altering the graphing parameters
  if (graph == "proportion") {
    x <- mat[,paste0(x.axis, ".fraction")]
    y <- mat[,paste0(y.axis, ".fraction")]
  } else if (graph == "count") {
    x <- mat[,x.axis]
    y <- mat[,y.axis] }
  if (dot.size != "total") {
    size <- mat[,dot.size]
  } else { size <- mat[,"sum"] }
  if (exportTable == TRUE) { return(mat) }
  
  plot <- ggplot(mat, aes(x=x, y = y, fill = class)) + 
                theme_classic() + 
                scale_fill_manual(values = .colorizer(palette,length(unique(mat$class)))) + 
                xlab(x.axis) + 
                ylab(y.axis) + 
                labs(size = "Total n")
  if (graph == "proportion") {
    plot <- plot + 
            geom_abline(slope = 1, intercept = 0, alpha = 0.4, lty=2)  + 
            scale_y_sqrt() + 
            scale_x_sqrt() 
  } else if (graph == "count") {
    plot <- plot + 
            ylim(0, max(x,y)) + 
            xlim(0, max(x,y)) 
  }
  plot <- plot + 
          geom_jitter(aes(size = size,), shape = 21, color = "black", stroke = 0.25)
  
  return(plot)  
}
