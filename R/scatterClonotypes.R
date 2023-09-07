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
#' scatterClonotype(combined, x.axis = "P17B", y.axis = "P17L",
#' graph = "proportion")
#' 
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param x.axis name of the list element to appear on the x.axis
#' @param y.axis name of the list element to appear on the y.axis
#' @param dot.size either total or the name of the list element to 
#' use for size of dots
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param graph graph either proportion or raw clonotype count
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @param seed the integer seed to set for the random variation of point coords.
#' 
#' @import ggplot2
#' 
#' @export
#' @return ggplot of the relative clonotype numbers

scatterClonotype <- function(df, 
                             cloneCall ="strict", 
                             x.axis = NULL, 
                             y.axis = NULL,
                             chain = "both",
                             dot.size = "total", 
                             split.by = NULL,
                             graph = "proportion", 
                             exportTable = FALSE,
                             palette = "inferno",
                             seed = NULL) {
  if (!is.null(seed)) {set.seed(seed)}
  
  df <- list.input.return(df, split.by)
  cloneCall <- theCall(cloneCall)
  axes <- which(names(df) %in% c(x.axis, y.axis, dot.size))
  if (chain != "both") {
    for (i in axes) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  x.df <- as.data.frame(table(df[[x.axis]][,cloneCall]))
  colnames(x.df)[2] <- x.axis
  y.df <- as.data.frame(table(df[[y.axis]][,cloneCall]))
  colnames(y.df)[2] <- y.axis
  combined.df <- merge(x.df, y.df, by = "Var1", all = TRUE)
  if (dot.size != "total") {
    if (dot.size %!in% colnames(combined.df)) {
      size.df <- as.data.frame(table(df[[dot.size]][,cloneCall]))
      colnames(size.df)[2] <- dot.size
      combined.df <- merge(combined.df, size.df, by = "Var1", all = TRUE) }
    combined.df[is.na(combined.df)] <- 0
    combined.df[,paste0("size", ".fraction")] <- combined.df[,dot.size]/sum(combined.df[,dot.size])
    labeling <- unique(c(x.axis, y.axis, dot.size))
  } else {
    combined.df[is.na(combined.df)] <- 0
    labeling <- unique(c(x.axis, y.axis)) }
  combined.df[,"class"] <- NA
  combined.df[,"sum"] <- rowSums(combined.df[,labeling])
  for (i in seq_along(labeling)) {
    if (length(labeling) > 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
    } else if (length(labeling) == 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
    }}
  combined.df[,"class"] <- ifelse(combined.df[,y.axis] >= 1 & combined.df[,x.axis] >= 1, paste0("dual.expanded"), combined.df[,"class"])
  combined.df[,paste0(x.axis, ".fraction")] <- combined.df[,x.axis]/sum(combined.df[,x.axis])
  combined.df[,paste0(y.axis, ".fraction")] <- combined.df[,y.axis]/sum(combined.df[,y.axis])
  if (graph == "proportion") {
    x <- combined.df[,paste0(x.axis, ".fraction")]
    y <- combined.df[,paste0(y.axis, ".fraction")]
  } else if (graph == "count") {
    x <- combined.df[,x.axis]
    y <- combined.df[,y.axis] }
  if (dot.size != "total") {
    size <- combined.df[,dot.size]
  } else { size <- combined.df[,"sum"] }
  if (exportTable == TRUE) { return(combined.df) }
  plot <- ggplot(combined.df, aes(x=x, y = y, fill = class)) + 
    theme_classic() + 
    scale_fill_manual(values = .colorizer(palette,length(unique(combined.df$class)))) + 
    xlab(x.axis) + ylab(y.axis) + labs(size = "Total n")
  if (graph == "proportion") {
    plot <- plot + geom_abline(slope = 1, intercept = 0, alpha = 0.4, lty=2)  + 
      scale_y_sqrt() + scale_x_sqrt() 
  } else if (graph == "count") {
    plot <- plot + ylim(0, max(x,y)) + xlim(0, max(x,y)) 
  }
  plot <- plot + geom_jitter(aes(size = size,), shape = 21, color = "black", stroke = 0.25)
  return(plot)  
}
