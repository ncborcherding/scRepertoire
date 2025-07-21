#' Calculate Startrac-based Diversity Indices 
#' 
#' @description This function utilizes the STARTRAC approach to calculate T cell 
#' diversity metrics based on the work of Zhang et al. (2018, Nature) 
#' [PMID: 30479382](https://pubmed.ncbi.nlm.nih.gov/30479382/). It can compute 
#' three distinct indices: clonal expansion (`expa`), cross-tissue migration 
#' (`migr`), and state transition (`tran`). 
#' 
#' @details
#' The function requires a `type` variable in the metadata, which specifies the
#' tissue origin or any other categorical variable for migration analysis.
#'
#' **Indices:**
#' \itemize{
#'   \item{\strong{expa (Clonal Expansion):}} Measures the extent of clonal 
#'         proliferation within a T cell cluster. It is calculated as 
#'         `1 - normalized Shannon entropy`. A higher value indicates greater 
#'         expansion of a few clones.
#'   \item{\strong{migr (Cross-Tissue Migration):}} Quantifies the movement of 
#'         clonal T cells across different tissues (as defined by the `type`
#'         parameter). It is based on the entropy of a clonotype's distribution 
#'         across tissues.
#'   \item{\strong{tran (State Transition):}} Measures the developmental 
#'         transition of clonal T cells between different functional clusters. 
#'         It is based on the entropy  of a clonotype's distribution across 
#'         clusters.
#' }
#'
#' **Pairwise Analysis:**
#' The `pairwise` parameter enables the calculation of migration or transition
#' between specific pairs of tissues or clusters, respectively.
#' \itemize{
#'   \item{For migration (`index = "migr"`), set `pairwise` to the `type` column 
#'         (e.g., `pairwise = "Type"`).
#'   \item{For transition (`index = "tran"`), set `pairwise` to `"majorCluster"`.}
#' }
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' scRep_example  <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
#' scRep_example$Type <- substring(scRep_example$orig.ident,4,4) 
#' 
#' # Calculate a single index (expansion)
#' StartracDiversity(scRep_example, 
#'                   type = "Type", 
#'                   group.by = "Patient",
#'                   index = "expa")
#'                   
#' # Calculate pairwise migration between tissue types
#' StartracDiversity(scRep_example, 
#'                   type = "Type", 
#'                   group.by = "Patient",
#'                   index = "migr",
#'                   pairwise = "Type")
#'
#' @param sc.data The single-cell object after [combineExpression()].
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC + nt). A custom column header can also be used.
#' @param chain The TCR/BCR chain to use. Use `both` to include both chains 
#' (e.g., TRA/TRB). Accepted values: `TRA`, `TRB`, `TRG`, `TRD`, `IGH`, `IGL` 
#' (for both light chains), `both`.
#' @param index A character vector specifying which indices to calculate. 
#' Options: "expa", "migr", "tran". Default is all three.
#' @param type The metadata variable that specifies tissue type for migration 
#' analysis.
#' @param group.by A column header in the metadata or lists to group the analysis 
#' by (e.g., "sample", "treatment"). If `NULL`, data will be analyzed as 
#' by list element or active identity in the case of single-cell objects.
#' @param pairwise The metadata column to be used for pairwise comparisons. 
#' Set to the `type` variable for pairwise migration or "majorCluster" for 
#' pairwise transition.
#' @param exportTable If `TRUE`, returns a data frame or matrix of the results 
#' instead of a plot.
#' @param palette Colors to use in visualization - input any 
#' [hcl.pals][grDevices::hcl.pals].
#' @param ... Additional arguments passed to the ggplot theme
#' @importFrom stats reshape
#' @export
#' @concept SC_Functions
#' @return A ggplot object visualizing STARTRAC diversity metrics or data.frame if
#'`exportTable = TRUE`.
#' @author Liangtao Zheng
StartracDiversity <- function(sc.data,
                              cloneCall = "strict", 
                              chain = "both",
                              index = c("expa", "migr", "tran"),
                              type = NULL,
                              group.by = NULL, 
                              pairwise = NULL,
                              exportTable = FALSE, 
                              palette = "inferno",
                              ...) {
  
  # Input validation
  index <- match.arg(index)
  if (!is.null(pairwise) && length(index) > 1) {
    stop("Pairwise analysis can only be performed for a single index ('migr' or 'tran').")
  }
  if (!is.null(pairwise) && !index %in% c("migr", "tran")) {
    stop("Pairwise analysis is only supported for 'migr' or 'tran' indices.")
  }
  
  
  # Prepare data
  df <- .grabMeta(sc.data)
  cloneCall <- .theCall(df, cloneCall)
  barcodes <- rownames(df)
  colnames(df)[ncol(df)] <- "majorCluster"
    
  if (is.null(group.by)) {
    if (!"orig.ident" %in% colnames(df)) {
      stop("Please select a group.by variable.")
    }
    group.by <- "orig.ident"
  }
  group.levels <- unique(df[,group.by])
  
  if (chain != "both") {
    df <- .offTheChain(df, chain, cloneCall)
  }

  # Process clonotypes
  df <- df %>%
    group_by(across(all_of(c(group.by, cloneCall)))) %>%
    dplyr::mutate(n = n()) %>%
    as.data.frame()
  
  rownames(df) <- barcodes
  remove.pos <- which(is.na(df[,cloneCall]) | df[,cloneCall] == "")
  if (length(remove.pos) > 0) {
    df <- df[-remove.pos,]
  }
  
  processed <- data.frame(
    Cell_Name = rownames(df),
    clone.id = df[,cloneCall],
    patient = df[,group.by],
    majorCluster = df[,"majorCluster"],
    loc = df[,type],
    stringsAsFactors = FALSE
  )
  processed[processed == "NA"] <- NA
  processed <- na.omit(processed)
  
  # Calculate indices
  mat.list <- lapply(group.levels, function(level) {
    subset_data <- processed[processed$patient == level,]
    if (!is.null(pairwise)) {
      comparison_col <- if (index == "migr") {
        "loc"
      } else { # index == "tran"
        "majorCluster"
      }
      .calculatePairwiseIndices(subset_data, index, comparison_col)
    } else {
      .calculateIndices(subset_data, index)
    }
  })
  
  mat <- bind_rows(mat.list, .id = "group")
  if (nrow(mat) == 0) {
    warning("No data available for calculation. Returning NULL.")
    return(NULL)
  }
  
  if(!is.null(pairwise)) {
    mat$variable <- index[1]
    mat <- mat[!is.nan(mat$value),]
  }
  
  if (exportTable) { 
    return(mat) 
  }
  # Plotting logic
  if (!is.null(pairwise)) {
    col_name <- colnames(mat)[grepl("comparison", colnames(mat))]
    plot <- ggplot(mat, aes(x = .data[[col_name]], y = .data$value)) +
      geom_boxplot(aes(fill = .data[[col_name]]), outlier.alpha = 0, na.rm = TRUE) +
      facet_wrap(~majorCluster) +
      labs(y = "Pairwise Index Score", x = "Comparison")
  } else {
    mat_melt <- reshape(mat,
                        varying = index,
                        v.names = "value",
                        timevar = "variable",
                        times = index,
                        direction = "long")
    values <- .alphanumericalSort(unique(mat_melt$majorCluster))
    mat_melt$majorCluster <- factor(mat_melt$majorCluster, levels = values)
    mat_melt$value <- as.numeric(mat_melt$value)
    col <- length(unique(mat_melt$majorCluster))
    
    plot <- ggplot(mat_melt, aes(x = majorCluster, y = .data[["value"]])) +
      geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0, na.rm = TRUE) +
      labs(y = "Index Score", 
           x = "Clusters") +
      theme(axis.title.x = element_blank())
  }
  
  plot <- plot + 
    .themeRepertoire(...) + 
    guides(fill = "none") +
    scale_fill_manual(values = .colorizer(palette, length(unique(mat$majorCluster))))
  
  if(length(index > 1)) {
    plot <- plot + facet_grid(variable ~ ., scales = "free_y") 
  }
  
  return(plot)
}


# Helper function for standard index calculation
.calculateIndices <- function(processed, indices) {
  if (nrow(processed) == 0) return(NULL)
  clonotype.dist.cluster <- table(processed[,c("clone.id", "majorCluster")])
  clonotype.dist.loc <- table(processed[,c("clone.id", "loc")])
  
  # Return NULL if no clusters are found
  if (ncol(clonotype.dist.cluster) == 0) return(NULL)
  
  calIndex.matrix <- data.frame(majorCluster = colnames(clonotype.dist.cluster))
  
  if ("expa" %in% indices) {
    entropy_val <- .mcolEntropy(clonotype.dist.cluster)
    entropy_max <- log2(colSums(clonotype.dist.cluster > 0))
    expa <- 1 - (entropy_val / entropy_max)
    calIndex.matrix$expa <- expa
  }
  
  # Check if there are clonotypes to process for migr/tran
  if (nrow(clonotype.dist.cluster) > 0) {
    clonotype.data <- data.frame(clone.id = rownames(clonotype.dist.cluster))
    weights.mtx <- sweep(clonotype.dist.cluster, 2, colSums(clonotype.dist.cluster), "/")
    
    if ("migr" %in% indices && "loc" %in% colnames(processed) && nrow(clonotype.dist.loc) > 0) {
      # Ensure clonotypes for migration calculation exist in the cluster distribution
      shared_clones_migr <- intersect(rownames(clonotype.dist.cluster), rownames(clonotype.dist.loc))
      if(length(shared_clones_migr) > 0) {
        clonotype.data$migr <- .mrowEntropy(clonotype.dist.loc[shared_clones_migr,,drop=FALSE])
        migr_matrix <- t(weights.mtx[shared_clones_migr, , drop=FALSE]) %*% as.matrix(clonotype.data$migr)
        calIndex.matrix$migr <- migr_matrix[,1]
      } else {
        calIndex.matrix$migr <- NA
      }
    }
    
    if ("tran" %in% indices) {
      clonotype.data$tran <- .mrowEntropy(clonotype.dist.cluster)
      tran_matrix <- t(weights.mtx) %*% as.matrix(clonotype.data$tran)
      calIndex.matrix$tran <- tran_matrix[,1]
    }
  } else {
    # If no clonotypes, set indices to NA
    if ("migr" %in% indices) calIndex.matrix$migr <- NA
    if ("tran" %in% indices) calIndex.matrix$tran <- NA
  }
  
  for (col in names(calIndex.matrix)) {
    if (is.numeric(calIndex.matrix[[col]])) {
      calIndex.matrix[[col]][is.nan(calIndex.matrix[[col]]) | is.infinite(calIndex.matrix[[col]])] <- NA
    }
  }
  
  return(calIndex.matrix)
}

# Helper function for pairwise index calculation
.calculatePairwiseIndices <- function(processed, index, pairwise_col) {
  if (nrow(processed) < 2) return(NULL)
  
  unique_items <- unique(processed[[pairwise_col]])
  if (length(unique_items) < 2) return(NULL)
  
  pairs <- combn(unique_items, 2, simplify = FALSE)
  
  pairwise_results <- lapply(pairs, function(p) {
    pair_data <- processed[processed[[pairwise_col]] %in% p,]
    
    if (index == "migr") {
      dist_table <- table(pair_data[,c("clone.id", "loc")])
      clonotype_dist_cluster <- table(pair_data[,c("clone.id", "majorCluster")])
    } else { # tran
      dist_table <- table(pair_data[,c("clone.id", "majorCluster")])
      clonotype_dist_cluster <- dist_table
    }
    
    if(nrow(dist_table) == 0) return(NULL)
    
    clonotype_data <- data.frame(clone.id = rownames(dist_table),
                                 value = .mrowEntropy(dist_table))
    
    weights_mtx <- sweep(clonotype_dist_cluster, 2, colSums(clonotype_dist_cluster), "/")
    
    # Ensure clone.ids match for matrix multiplication
    shared_clones <- intersect(rownames(weights_mtx), clonotype_data$clone.id)
    if (length(shared_clones) == 0) return(NULL)
    
    weights_mtx_filtered <- weights_mtx[shared_clones,, drop=FALSE]
    clonotype_data_filtered <- clonotype_data[clonotype_data$clone.id %in% shared_clones,]
    
    result_matrix <- t(weights_mtx_filtered) %*% as.matrix(clonotype_data_filtered$value)
    
    res <- data.frame(majorCluster = rownames(result_matrix), value = result_matrix[,1])
    res$comparison <- paste(p, collapse = " vs ")
    return(res)
  })
  
  bind_rows(pairwise_results)
}

# Entropy of each row of the input matrix
.mrowEntropy <- function(x) {
  freqs <- sweep(x, 1, rowSums(x), "/")
  H <- -rowSums(ifelse(freqs > 0, freqs * log2(freqs), 0))
  return(H)
}

# Entropy of each column of the input matrix
.mcolEntropy <- function(x) {
  freqs <- sweep(x, 2, colSums(x), "/")
  H <- -colSums(ifelse(freqs > 0, freqs * log2(freqs), 0))
  return(H)
}

