#' Examining the VDJ gene usage across clones
#'
#' This function the proportion, percent, or counts of V, D, or J genes used by 
#' grouping variables. This function only quantifies single gene loci for 
#' indicated **chain**. For examining VJ or other gene parings, please 
#' see [percentVJ()] or [vizGenes()].
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' percentGenes(combined, 
#'              chain = "TRB", 
#'              gene = "Vgene")
#' 
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param chain indicate a specific chain should be used - 
#' e.g. "TRA", "TRG", "IGH", "IGL", etc
#' @param gene "V", "D" or "J"
#' @param group.by The variable to use for grouping
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order
#' @param summary.fun Character string choosing the summary statistic - 
#' `"percent"` (default), `"proportion"`, or `"count"`
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#' 
#' @importFrom immApex calculateGeneUsage
#' 
#' @export
#' 
#' @concept Summarize_Repertoire
#' 
#' @return A ggplot object displaying a heatmap of gene percentages.
#' If `exportTable = TRUE`, a matrix of the raw data is returned.
percentGenes <- function(input.data,
                         chain = "TRB",
                         gene = "Vgene", 
                         group.by = NULL, 
                         order.by = NULL,
                         exportTable = FALSE,
                         summary.fun = c("percent", "proportion", "count"),
                         palette = "inferno") {
  
  sco <- .is.seurat.or.se.object(input.data)
  summary.fun <- match.arg(summary.fun)
  input.data <- .dataWrangle(input.data, group.by, "CTgene", chain)
  if(!is.null(group.by) & !sco) {
    input.data <- .groupList(input.data, group.by)
  }
  
  #Parsing gene input
  gene_map <- c(Vgene = "V", V = "V", v = "V", "v.gene" = "V",
                Jgene = "J", j = "J", J = "J", "j.gene" = "J",
                Dgene = "D", d = "D", D = "D", "D.gene" = "D")
  gene_type <- gene_map[gene]
  
  if (is.na(gene_type)) {
    stop("Invalid gene argument. Please use 'V', 'D', or 'J'.")
  }
  
  if (gene_type == "D" && !(chain %in% c("TRB", "TRD", "IGH"))) {
    stop(paste0("There is no D locus for ", chain, "."))
  }
  
  #Getting genes counts across list
  gene.counts <- lapply(input.data, function(x) {
    tmp <- unlist(str_split(x[,"CTgene"], ";"))
    tmp <- str_split(tmp, "[.]", simplify = TRUE)
    tmp <- tmp[!is.na(tmp[,1]),]
    colnames(tmp) <- c("V", "D", "J", "C")
    counts <- immApex::calculateGeneUsage(as.data.frame(tmp), 
                                          loci = gene_type, 
                                          summary = summary.fun)
  })
  
  unique.genes <- unique(unlist(lapply(gene.counts, names)))
  mat <- matrix(0, ncol = length(unique.genes), nrow = length(input.data))
  colnames(mat) <- unique.genes
  rownames(mat) <- names(input.data)
  
  # Populating matrix
  for (i in seq_along(gene.counts)) {
    sample_counts <- gene.counts[[i]]
    col_indices <- match(names(sample_counts), colnames(mat))
    mat[i, na.omit(col_indices)] <- sample_counts[!is.na(col_indices)]
  }
  
  if (exportTable == TRUE) { 
    return(mat) 
  }
  
  # Getting mat into a ggplot-compliant form
  mat_melt <- expand.grid(Var1 = rownames(mat),
                          Var2 = colnames(mat), 
                          KEEP.OUT.ATTRS = FALSE)
  mat_melt$value <- as.vector(mat)
  if(!is.null(order.by)) {
    mat_melt <- .orderingFunction(vector = order.by,
                                  group.by = "Var1", 
                                  mat_melt)
  }
  
  plot <- ggplot(mat_melt, aes(y=Var1, x = Var2, fill=value)) +
    geom_tile(lwd= 0.25, color = "black") +
    scale_fill_gradientn(name = .toCapitilize(summary.fun), 
                         colors = .colorizer(palette,21)) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          axis.title = element_blank())
  return(plot)
}
