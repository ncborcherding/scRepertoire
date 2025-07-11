#' Summarizes and Visualizes Gene Usage
#'
#' This function quantifies and visualizes the usage of V, D, or J genes,
#' or gene pairings across T or B cell clones. 
#' 
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#'                                     
#' # Visualize single gene (TRBV) usage as a heatmap, grouped by sample
#' percentGeneUsage(combined.TCR,
#'                  genes = "TRBV",
#'                  group.by = "sample",
#'                  plot.type = "heatmap",
#'                  summary.fun = "percent")
#'
#' # Visualize single gene (TRBV) usage as a barplot, grouped by sample
#' percentGeneUsage(combined.TCR,
#'                  genes = "TRBV",
#'                  group.by = "sample",
#'                  plot.type = "barplot",
#'                  summary.fun = "count")
#'
#' # Visualize paired gene (TRBV-TRBJ) usage as a heatmap
#' percentGeneUsage(combined.TCR,
#'                  genes = c("TRBV", "TRBJ"),
#'                  group.by = "sample",
#'                  plot.type = "heatmap",
#'                  summary.fun = "proportion")
#'
#' # Export the raw data table for single gene usage
#' trbv_usage_table <- percentGeneUsage(combined.TCR,
#'                                      genes = "TRBV",
#'                                      group.by = "sample",
#'                                      exportTable = TRUE,
#'                                      summary.fun = "count")
#'
#' # Export the raw data table for paired gene usage
#' trbv_trbj_usage_table <- percentGeneUsage(combined.TCR,
#'                                           genes = c("TRBV", "TRBJ"),
#'                                           group.by = "sample",
#'                                           exportTable = TRUE,
#'                                           summary.fun = "percent")
#' 
#'
#' @param input.data The product of [combineTCR()], [combineBCR()],
#'or [combineExpression()].
#' @param genes A character vector specifying the gene loci to analyze.
#' Can be a single gene e.g., "TRBV" or "IGHJ" or a pair for genes analysis 
#' (e.g., c("TRBV", "TRAV"), or "TRBV", "TRBJ").
#' @param group.by The variable to use for grouping (e.g., "sample", "orig.ident").
#' @param order.by A vector of specific plotting order or "alphanumeric"
#' to plot groups in order.
#' @param summary.fun Character string choosing the summary statistic -
#' `"percent"` (default), `"proportion"`, or `"count"`.
#' @param plot.type The type of plot to return: `"heatmap"` (default for paired loci,
#' also available for single loci), or `"barplot"` (for single loci).
#' @param exportTable Logical. If `TRUE`, returns the data frame or matrix used
#' for forming the graph instead of the plot.
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#'
#' @importFrom immApex calculateGeneUsage getIR
#' @importFrom graphics plot
#' @importFrom grDevices hcl.pals
#' @importFrom utils stack
#' @importFrom dplyr group_by summarise ungroup mutate 
#'
#' @export
#'
#' @concept Summarize_Repertoire
#' @concept Visualizing_Clones
#'
#' @return A ggplot object displaying a heatmap or bar plot of gene usage.
#' If `exportTable = TRUE`, a matrix or data frame of the raw data is returned.
percentGeneUsage <- function(input.data,
                             chain = "TRB",
                             genes = "TRBV",
                             group.by = NULL,
                             order.by = NULL,
                             summary.fun = c("percent", "proportion", "count"),
                             plot.type = "heatmap", 
                             exportTable = FALSE,
                             palette = "inferno") {
  
  sco <- .is.seurat.or.se.object(input.data)
  
  summary.fun <- match.arg(summary.fun)
  if (!is.character(genes) || !(length(genes) %in% 1:2)) {
    stop("Parameter 'genes' must be a character vector of length 1 or 2, specifying gene loci.")
  }
  
  chains_to_extract <- unique(substr(genes, 1,3))
  if (length(chains_to_extract) > 2) {
    stop("Analysis limited to a maximum of two chains. Please refine 'genes' input.")
  }
  
  gene_segment_map <- c("V" = "v", "D" = "d", "J" = "j", "C" = "c")
  column.headers <- character(length(genes))
  for (i in seq_along(genes)) {
    gene_type_char <- substr(genes[i], 4, 4)
    if (gene_type_char %in% names(gene_segment_map)) {
      column.headers[i] <- gene_segment_map[gene_type_char]
    } else {
      stop(paste0("Unsupported gene segment format: '", genes[i], "'. Expected format like 'TRBV', 'IGHJ', etc."))
    }
  }
  
  ir_data_list <- lapply(chains_to_extract, function(ch) {
    dfs <- immApex::getIR(input.data, 
                          chains = ch, 
                          sequence.type = "aa", 
                          group.by = group.by)
    if(inherits(dfs, "list")) {
      dfs <- do.call(rbind, dfs) 
    } 
    return(dfs)
  })
  
  # Handling NULL group.by
  if(is.null(group.by)) {
    if(sco) {
      group.by <- "ident"
      groupings <- .bound.input.return(input.data)[,group.by, drop = FALSE]
    } else {
      group.by <- "element.names" 
      groupings <- .bound.input.return(input.data)[,c("barcode", group.by), drop = FALSE]
      rownames(groupings) <- groupings[["barcode"]]
      groupings <- groupings[,-which(colnames(groupings) == "barcode"), drop = FALSE]
      
    }
    ir_data_list <- lapply(ir_data_list, function(x) {
      merge(x, groupings, by.x = "barcode", by.y = 0)
    })
  }
  gene.data <- ir_data_list[[1]]

  
  # If second chain is being used
  if (length(chains_to_extract) == 2 && length(genes) == 2) {
    second_chain_data <- ir_data_list[[2]]
    common_cols <- intersect(names(gene.data), names(second_chain_data))
    cols_to_rename <- setdiff(names(second_chain_data), c(group.by, "barcode"))
    colnames(second_chain_data)[match(cols_to_rename, colnames(second_chain_data))] <- paste0(cols_to_rename, ".2")
    gene.data <- merge(gene.data, second_chain_data, by = c("barcode", group.by), all = TRUE)
    column.headers[2] <- paste0(column.headers[2], ".2")
  }
  
  # Final separation before calculation
  grouped_data <- split(gene.data, gene.data[[group.by]])
  calculated_usage <- lapply(grouped_data, function(x) {
      df_for_calc <- x
      immApex::calculateGeneUsage(
        input.data = df_for_calc,
        loci = column.headers,
        summary = summary.fun
      )
    })
  
  mat_melt_list <- lapply(names(calculated_usage), function(group_name) {
    item <- calculated_usage[[group_name]]
    if (is.atomic(item) && !is.null(names(item)) && length(dim(item)) < 2) {
      # Single locus result (vector)
      data.frame(
        Var1 = names(item),
        Weight = as.vector(unname(item)),
        Group = group_name,
        stringsAsFactors = FALSE
      )
    } else if (is.matrix(item)) {
      # Paired loci result 
      df_melted <- as.data.frame(as.table(item), stringsAsFactors = FALSE)
      names(df_melted) <- c("Var1", "Var2", "Weight")
      data.frame(
        Var1 = df_melted$Var1,
        Var2 = df_melted$Var2,
        Weight = df_melted$Weight,
        Group = group_name,
        stringsAsFactors = FALSE
      )
    }
  })
  
  mat_melt <- do.call(rbind, mat_melt_list)
  if (exportTable) {
    if (length(genes) == 1) {
      # Single gene usage: rows are genes, columns are groups
      output_matrix <- tapply(mat_melt$Weight, list(mat_melt$Var1, mat_melt$Group), sum)
    } else {
      # Paired gene usage: rows are gene pairs, columns are groups
      row_variable <- paste(mat_melt$Var1, mat_melt$Var2, sep = "_")
      output_matrix <- tapply(mat_melt$Weight, list(row_variable, mat_melt$Group), sum)
    }
    output_matrix[is.na(output_matrix)] <- 0
    return(output_matrix)
  }
  
  col.lab <- .toCapitilize(summary.fun)
  
  if (!is.null(order.by)) {
      mat_melt <- .orderingFunction(vector = order.by, group.by = "Group", mat_melt)
  }
  
  if (length(genes) == 1) {
    if (plot.type == "barplot") {
      plot <- ggplot(mat_melt, aes(x = Var1, y = Weight)) +
        geom_bar(stat = "identity", color = "black", lwd = 0.25) +
        theme_classic() +
        labs(y = col.lab) + 
        theme(axis.title.x = element_blank(),,
              axis.ticks.x = element_blank(),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=rel(0.8)))
      
      if (length(unique(mat_melt$Group)) > 1) {
        plot <- plot + facet_grid(Group ~ .)
      }
      
    } else { 
      plot <- ggplot(mat_melt, aes(y = Var1, x = Group, fill = Weight)) +
        geom_tile(lwd = 0.1, color = "black") +
        scale_fill_gradientn(colors = .colorizer(palette, 21)) +
        labs(fill = col.lab) + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
              axis.title = element_blank())
    }
  } else { # Paired genes (heatmap only)
    plot <- ggplot(mat_melt, aes(y = Var1, x = Var2, fill = round(Weight, 2))) +
      geom_tile(lwd = 0.1, color = "black") +
      scale_fill_gradientn(colors = .colorizer(palette, 21)) +
      labs(fill = col.lab) + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title = element_blank())
    if (length(unique(mat_melt$Group)) > 1) {
      plot <- plot + facet_wrap(~Group)
    }
  }
  return(plot)
}

#' @rdname percentGeneUsage
#' @examples
#' # Visualize paired gene (TRBV-TRBJ) usage as a heatmap
#' vizGenes(combined.TCR,
#'          x.axis = "TRBV",
#'          y.axis = "TRBJ",
#'          group.by = "sample",
#'          summary.fun = "count")
#'
#' # Visualize cross-chain gene pairing (TRBV-TRAV)
#' vizGenes(combined.TCR,
#'          x.axis = "TRBV",
#'          y.axis = "TRAV",
#'          group.by = "sample",
#'          summary.fun = "percent")
#' 
#' @export
vizGenes <- function(input.data,
                     x.axis = "TRBV",
                     y.axis = NULL, 
                     group.by = NULL,
                     plot = "heatmap",
                     order.by = NULL, 
                     summary.fun = c("percent", "proportion", "count"),
                     exportTable = FALSE,
                     palette = "inferno") {
  
  summary.fun <- match.arg(summary.fun)
  genes_param <- x.axis
  plot_type_param <- plot
  group_by_param <- group.by
  order_by_param <- order.by 
  
  # Handle y.axis for paired gene analysis
  if (!is.null(y.axis) && grepl("TRA|TRB|TRG|TRD|IGH|IGL|IGK", y.axis)) {
    genes_param <- c(x.axis, y.axis)
    plot_type_param <- "heatmap" 
  } else if (!is.null(y.axis)) {
    # If y.axis is a categorical variable, use it as group.by
    group_by_param <- y.axis
  } else if(is.null(y.axis)) {
    y.axis <- group_by_param
  }
  
  percentGeneUsage(
    input.data = input.data,
    genes = genes_param,
    group.by = group_by_param,
    order.by = order_by_param,
    summary.fun = summary.fun, 
    plot.type = plot_type_param,
    exportTable = exportTable,
    palette = palette
  )
}

#' @rdname percentGeneUsage
#' @examples
#' 
#' # Quantify and visualize TRA V-gene usage as a heatmap
#' percentGenes(combined.TCR,
#'              chain = "TRA",
#'              gene = "Vgene",
#'              group.by = "sample",
#'              summary.fun = "percent")
#'
#' # Quantify TRA J-gene usage and export the table
#' ighj_usage_table <- percentGenes(combined.TCR,
#'                                  chain = "TRA",
#'                                  gene = "Jgene",
#'                                  group.by = "sample",
#'                                  exportTable = TRUE,
#'                                  summary.fun = "count")
#' 
#' @export
percentGenes <- function(input.data,
                         chain = "TRB",
                         gene = "Vgene",
                         group.by = NULL,
                         order.by = NULL,
                         exportTable = FALSE,
                         summary.fun = c("percent", "proportion", "count"),
                         palette = "inferno") {
  
  summary.fun <- match.arg(summary.fun)
  
  
  gene_segment_prefix <- switch(chain,
                                "TRA" = "TRA", "TRB" = "TRB", "TRG" = "TRG", "TRD" = "TRD",
                                "IGH" = "IGH", "IGL" = "IGL", "IGK" = "IGK",
                                stop("Invalid chain argument.")
  )
  
  gene_type_suffix <- switch(tolower(substr(gene,1,1)),
                             "v" = "V", "d" = "D", "j" = "J",
                             stop("Invalid gene argument. Please use 'V', 'D', or 'J'.")
  )
  
  # Construct the 'genes' argument for percentGeneUsage
  genes_param <- paste0(gene_segment_prefix, gene_type_suffix)
  
  percentGeneUsage(
    input.data = input.data,
    genes = genes_param,
    group.by = group.by,
    order.by = order.by,
    summary.fun = summary.fun,
    plot.type = "heatmap", 
    exportTable = exportTable,
    palette = palette
  )
}

#' @rdname percentGeneUsage
#' @examples

#' # Quantify and visualize TRB V-J gene pairings as a heatmap
#' percentVJ(combined.TCR,
#'           chain = "TRB",
#'           group.by = "sample",
#'           summary.fun = "percent")
#'
#' # 2. Quantify TRA V-J gene pairings and export the table
#' trav_traj_table <- percentVJ(combined.TCR,
#'                              chain = "TRA",
#'                              group.by = "sample",
#'                              exportTable = TRUE,
#'                              summary.fun = "proportion")
#' 
#' @export
percentVJ <- function(input.data,
                      chain = "TRB",
                      group.by = NULL,
                      order.by = NULL,
                      summary.fun = c("percent", "proportion", "count"),
                      exportTable = FALSE,
                      palette = "inferno") {
  
  summary.fun <- match.arg(summary.fun)
  genes_param <- switch(chain,
                        "TRA" = c("TRAV", "TRAJ"),
                        "TRB" = c("TRBV", "TRBJ"),
                        "TRG" = c("TRGV", "TRGJ"),
                        "TRD" = c("TRDV", "TRDJ"), 
                        "IGH" = c("IGHV", "IGHJ"),
                        "IGL" = c("IGLV", "IGLJ"),
                        "IGK" = c("IGKV", "IGKJ"),
                        stop("Invalid chain argument.")
  )
  
  percentGeneUsage(
    input.data = input.data,
    genes = genes_param,
    group.by = group.by,
    order.by = order.by,
    summary.fun = summary.fun, 
    plot.type = "heatmap",
    exportTable = exportTable,
    palette = palette
  )
}
