#' Calculate Clonal Diversity
#'
#' This function calculates a specified diversity metric for samples or groups
#' within a dataset. To control for variations in library size, the function
#' can perform bootstrapping with downsampling. It resamples each group to the
#' size of the smallest group and calculates the diversity metric across
#' multiple iterations, returning the mean value.
#' 
#' @details
#' The function operates by first splitting the dataset by the specified `group.by`
#' variable.
#'
#' **Downsampling and Bootstrapping:**
#' To make a fair comparison between groups of different sizes, diversity metrics
#' often require normalization. This function implements this by downsampling.
#' 1.  It determines the number of clones in the smallest group.
#' 2.  For each group, it performs `n.boots` iterations (default = 100).
#' 3.  In each iteration, it randomly samples the clones (with replacement) down to
#'     the size of the smallest group.
#' 4.  It calculates the selected diversity metric on this downsampled set.
#' 5.  The final reported diversity value is the mean of the results from all
#'     bootstrap iterations.
#'
#' This process can be disabled by setting `skip.boots = TRUE`.
#' 
#' Available Diversity Metrics (metric):
#' The function uses a registry of metrics imported from the immApex package.
#' You can select one of the following:
#' \itemize{
#'   \item{\code{"shannon"}: Shannon's Entropy. See \code{\link[immApex]{shannon_entropy}}.}
#'   \item{\code{"inv.simpson"}: Inverse Simpson Index. See \code{\link[immApex]{inv_simpson}}.}
#'   \item{\code{"gini.simpson"}: Gini-Simpson Index. See \code{\link[immApex]{gini_simpson}}.}
#'   \item{\code{"norm.entropy"}: Normalized Shannon Entropy. See \code{\link[immApex]{norm_entropy}}.}
#'   \item{\code{"pielou"}: Pielou's Evenness (same as norm.entropy). See \code{\link[immApex]{pielou_evenness}}.}
#'   \item{\code{"ace"}: Abundance-based Coverage Estimator. See \code{\link[immApex]{ace_richness}}.}
#'   \item{\code{"chao1"}: Chao1 Richness Estimator. See \code{\link[immApex]{chao1_richness}}.}
#'   \item{\code{"gini"}: Gini Coefficient for inequality. See \code{\link[immApex]{gini_coef}}.}
#'   \item{\code{"d50"}: The number of top clones making up 50% of the library. See \code{\link[immApex]{d50_dom}}.}
#'   \item{\code{"hill0"}, \code{"hill1"}, \code{"hill2"}: Hill numbers of order 0, 1, and 2. See \code{\link[immApex]{hill_q}}.}
#' }
#'
#' @examples
#' # Making combined contig data
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B","P19L", "P20B", "P20L"))
#'
#' # Calculate Shannon diversity, grouped by sample
#' clonalDiversity(combined, 
#'                 cloneCall = "gene", 
#'                 metric = "shannon")
#'
#' # Calculate Inverse Simpson without bootstrapping
#' clonalDiversity(combined, 
#'                 cloneCall = "aa", 
#'                 metric = "inv.simpson", 
#'                 skip.boots = TRUE)
#'
#' @param input.data The product of [combineTCR()],
#' [combineBCR()], or [combineExpression()].
#' @param cloneCall How to call the clone - VDJC gene (**gene**),
#' CDR3 nucleotide (**nt**), CDR3 amino acid (**aa**),
#' VDJC gene + CDR3 nucleotide (**strict**), or a custom variable
#' in the data.
#' @param metric The diversity metric to calculate. Must be a single string from
#' the list of available metrics (see Details).
#' @param chain Indicate if both or a specific chain should be used -
#' e.g. "both", "TRA", "TRG", "IGH", "IGL".
#' @param group.by Variable in the metadata to group samples for calculation.
#' @param order.by A vector of specific plotting order for the `group.by` variable,
#' or "alphanumeric" to plot groups in that order.
#' @param x.axis An additional metadata variable to group samples along the x-axis.
#' @param exportTable If TRUE, returns a data frame of the results instead of a plot.
#' @param return.boots If TRUE, returns all bootstrap values instead of the mean.
#' Automatically sets `exportTable = TRUE`.
#' @param skip.boots If TRUE, disables downsampling and bootstrapping. The metric
#' will be calculated on the full dataset for each group.
#' @param n.boots The number of bootstrap iterations to perform (default is 100).
#' @param palette Colors to use in visualization - input any
#' [hcl.pals][grDevices::hcl.pals].
#'
#' @import ggplot2
#' @export
#' @concept Visualizing_Clones
#' @return A ggplot object visualizing the diversity metric, or a data.frame if
#' `exportTable = TRUE`.
#' @author Andrew Malone, Nick Borcherding, Nathan Vanderkraan
clonalDiversity <- function(input.data,
                            cloneCall = "strict",
                            metric = "shannon",
                            chain = "both",
                            group.by = NULL,
                            order.by = NULL,
                            x.axis = NULL,
                            exportTable = FALSE,
                            palette = "inferno",
                            n.boots = 100,
                            return.boots = FALSE,
                            skip.boots = FALSE) {
  
  #Argument and Data Validation 
  if (length(metric) != 1 || !is.character(metric)) {
    stop("`metric` must be a single string.")
  }
  if (!metric %in% names(.div.registry)) {
    stop("Invalid `metric`. Please choose from: ", paste(names(.div.registry), collapse = ", "))
  }
  if(return.boots) {
    exportTable <- TRUE
  }
  
  # Data Wrangling 
  sco <- .is.seurat.or.se.object(input.data)
  cloneCall <- .theCall(input.data, cloneCall, check.df = FALSE)
  input.data <- .dataWrangle(input.data, group.by, cloneCall, chain)
  
  if(!is.null(group.by) && !sco) {
    grouping <- c(group.by, x.axis)
    input.data <- .groupList(input.data, grouping)
  }
  
  # Core Calculation 
  div_func <- .div.registry[[metric]]
  min_size <- min(sapply(input.data, function(x) nrow(x)))
  
  # Efficiently process each group using lapply
  results_list <- lapply(names(input.data), function(group_name) {
    sub_data <- input.data[[group_name]]
    clones <- sub_data[[cloneCall]]
    
    if (skip.boots) {
      diversity_scores <- div_func(table(clones))
    } else {
      # Use replicate for efficient bootstrapping
      diversity_scores <- replicate(n.boots, {
        resampled_clones <- sample(clones, size = min_size, replace = TRUE)
        div_func(table(resampled_clones))
      })
    }
    
    if (return.boots) {
      data.frame(
        group = group_name,
        value = diversity_scores,
        metric = metric
      )
    } else {
      data.frame(
        group = group_name,
        value = mean(diversity_scores, na.rm = TRUE),
        metric = metric
      )
    }
  })
  
  # Combine results into a single data frame 
  output_df <- do.call(rbind, results_list)
  colnames(output_df)[colnames(output_df) == "group"] <- group.by %||% "Group"
  
  # Prepare for plotting or export 
  if (!is.null(x.axis)) {
    x.variable <- lapply(input.data, function(x) {
      unique(x[,x.axis])[1]
    })
    x.variable <- do.call(rbind, x.variable)
    x.variable <- as.data.frame(x.variable)
    colnames(x.variable) <- x.axis
    output_df <- merge(output_df, x.variable, by.x = group.by %||% "Group", by.y = "row.names")
    output_df[,1] <- gsub("\\..*", "", output_df[,1]) #removing the x.axis grouping
  } else {
    x.axis <- "x.axis"
    output_df[,x.axis] <- 1
  }
  
  if (exportTable) {
    return(output_df)
  }
  
  # Arranging order if using order.by
  if(!is.null(order.by)) {
    output_df <- .orderingFunction(vector = order.by,
                                   group.by = names(output_df)[1],
                                   mat_melt = output_df)
  }
  
  # Plotting 
  metric.name <- gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", metric, perl=TRUE)
  plot <- ggplot(output_df, aes(x = .data[[x.axis]], y = as.numeric(value))) +
    geom_boxplot(outlier.alpha = 0) +
    geom_jitter(aes(fill = .data[[group.by %||% "Group"]]),
                size = 3,
                shape = 21,
                stroke = 0.25,
                color = "black") +
    labs(fill = "Group", y = paste(metric.name, "Index Score")) +
    scale_fill_manual(values = .colorizer(palette, length(unique(output_df[[group.by %||% "Group"]])))) +
    theme_classic() +
    theme(axis.title.x = element_blank())
  
  if (x.axis == "x.axis") {
    plot <- plot + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  }
  
  return(plot)
}

# Helper for null default
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' @importFrom immApex ace_richness chao1_richness gini_coef d50_dom 
#' shannon_entropy inv_simpson gini_simpson norm_entropy pielou_evenness hill_q
.div.registry <- list(
  ace          = ace_richness,
  chao1        = chao1_richness,
  gini         = gini_coef,
  d50          = d50_dom,
  shannon      = shannon_entropy,
  inv.simpson  = inv_simpson,
  gini.simpson = gini_simpson,
  norm.entropy = norm_entropy,
  pielou       = pielou_evenness,
  hill0        = hill_q(0),   # richness
  hill1        = hill_q(1),   # exp(H)
  hill2        = hill_q(2)    # 1/Simpson
)
