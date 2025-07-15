#' Compare Clonal Abundance Across Variables
#'
#' This function visualizes the relative abundance of specific clones across 
#' different samples or groups. It is useful for tracking how the proportions 
#' of top clones change between conditions. The output can be an alluvial plot 
#' to trace clonal dynamics or an area plot to show compositional changes.
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B","P19L", "P20B", "P20L"))
#' 
#' # Using clonalCompares()
#' clonalCompare(combined,
#'               top.clones = 5,
#'               samples = c("P17B", "P17L"),
#'               cloneCall="aa")
#'
#' @param input.data The product of \code{\link{combineTCR}},
#' \code{\link{combineBCR}}, or \code{\link{combineExpression}}.
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param samples The specific samples to isolate for visualization.
#' @param clones The specific clonal sequences of interest
#' @param top.clones The top number of clonal sequences per group.
#' (e.g., top.clones = 5)
#' @param highlight.clones Clonal sequences to highlight, if present,
#' all other clones returned will be grey
#' @param relabel.clones Simplify the legend of the graph by returning
#' clones that are numerically indexed
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed 
#' by list element or active identity in the case of single-cell objects.
#' @param order.by A character vector defining the desired order of elements 
#' of the `group.by` variable. Alternatively, use `alphanumeric` to sort groups 
#' automatically.
#' @param graph The type of plot to generate. Accepted values are `alluvial` 
#' (default) or `area`
#' @param proportion If `TRUE` (default), the y-axis will represent the 
#' proportional abundance of clones. If `FALSE`, the y-axis will represent 
#' raw clone counts.`
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any
#' \link[grDevices]{hcl.pals}
#' @param ... Additional arguments passed to the ggplot theme
#' 
#' @export
#' @importFrom dplyr slice_max
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing proportions of clones by groupings, or a
#' data.frame if `exportTable = TRUE`.
clonalCompare <- function(input.data,
                          cloneCall = "strict",
                          chain = "both",
                          samples = NULL,
                          clones = NULL,
                          top.clones = NULL,
                          highlight.clones = NULL,
                          relabel.clones = FALSE,
                          group.by = NULL,
                          order.by = NULL,
                          graph = "alluvial",
                          proportion = TRUE,
                          exportTable = FALSE,
                          palette = "inferno",
                          ...) {

  #Tie goes to indicated clones over top clones
  if(!is.null(top.clones) && !is.null(clones)) {
    top.clones <- NULL
  }
  input.data <- .dataWrangle(input.data,
                              group.by,
                              .theCall(input.data, cloneCall, check.df = FALSE),
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)

  sco <- .is.seurat.or.se.object(input.data)
  if(!is.null(group.by) && !sco) {
    input.data <- .groupList(input.data, group.by)
  }

  compareColname <- ifelse(proportion, "Proportion", "Count")
  normalizer <- ifelse(proportion, sum, length)

  Con.df <- input.data %>%
    purrr::imap(function(df, columnNames) {
      tbl <- as.data.frame(table(df[, cloneCall]))
      if (proportion) {
        tbl[, 2] <- tbl[, 2] / normalizer(tbl[, 2])
      }
      colnames(tbl) <- c("clones", compareColname)
      tbl$Sample <- columnNames
      tbl
    }) %>%
    dplyr::bind_rows()

  #Filtering steps
  if (!is.null(samples)) {
    Con.df <- Con.df[Con.df$Sample %in% samples,]
  }
  if (!is.null(clones)) {
    Con.df <- Con.df[Con.df$clones %in% clones,]
  } else if (!is.null(top.clones)) {
    top <- Con.df %>%
      group_by(Sample) %>%
      slice_max(
        n = top.clones,
        order_by = !!sym(compareColname),
        with_ties = FALSE
      )
    Con.df <- Con.df[Con.df$clones %in% top$clones,]
  }
  if (nrow(Con.df) < length(unique(Con.df$Sample)) || nrow(Con.df) == 0) {
    stop("Please reasses the filtering strategies here, there are not
            enough clones to examine.")
  }
  #Clones relabeling
  clones.returned <- as.vector(unique(Con.df[order(Con.df[, compareColname], decreasing = TRUE),"clones"]))
  if (relabel.clones) {
    new.clones <- paste0("Clone: ", seq_len(length(clones.returned)))
    names(new.clones) <- clones.returned
    #Isolated new clone names for highlight purposes
    if(!is.null(highlight.clones)) {
      highlight.clones <- unname(new.clones[which(names(new.clones) %in% highlight.clones)])
    }
    Con.df[,"original.clones"] <- Con.df[, "clones"]
    Con.df[,"clones"] <- new.clones[as.vector(Con.df[, "clones"])]
    Con.df[,"clones"] <- factor(Con.df[, "clones"],
                                levels = .alphanumericalSort(unique(Con.df[, "clones"])))
    clones.returned <- as.vector(unique(Con.df[, "clones"]))
  }

  if(!is.null(order.by)) {
    Con.df <- .orderingFunction(vector = order.by,
                                group.by = "Sample",
                                data.frame = Con.df)
  }

  if (exportTable) {
    return(Con.df)
  }

  #Plotting Functions
  plot <- ggplot(Con.df, aes(x = .data[["Sample"]],
                             fill = clones,
                             group = clones,
                             stratum = clones,
                             alluvium = clones,
                             y = !!sym(compareColname),
                             label = clones)) +
    .themeRepertoire(...) +
    theme(axis.title.x = element_blank(),
          legend.text=element_text(size=rel(0.5)),
          legend.key.size = unit(0.5,"line"))
  if (graph == "alluvial") {
    plot <- plot +  geom_stratum() + geom_flow(stat = "alluvium")
  } else if (graph == "area") {
    plot <- plot +
      geom_area(aes(group = clones), color = "black")
  }

  #Highlighting specific clones
  if (!is.null(highlight.clones)) {
    clone.colors <- rep("grey", length(clones.returned))
    pos <- which(clones.returned %in% highlight.clones)
    clone.colors[pos] <- .colorizer(palette, length(pos))
    names(clone.colors) <- clones.returned
    plot <- plot + scale_fill_manual(values = clone.colors)
  } else {
    plot <- plot + scale_fill_manual(values = .colorizer(palette, length(unique(Con.df[,"clones"]))))
  }
  return(plot)
}
