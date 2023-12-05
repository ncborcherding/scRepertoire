#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of 
#' various overlap metrics for clonotypes.
#' 
#' @details
#' The formulas for the indices are as follows:
#' 
#' \strong{Overlap Coefficient:}
#' \deqn{overlap = \frac{\sum \min(a_i, b_i)}{\min(\sum a_i, \sum b_i)}}
#' 
#' \strong{Raw Count Overlap:}
#' \deqn{raw = \sum \min(a_i, b_i)}
#' 
#' \strong{Morisita Index:}
#' \deqn{morisita = \frac{\sum a_i b_i}{(\sum a_i)(\sum b_i)}}
#' 
#' \strong{Jaccard Index:}
#' \deqn{jaccard = \frac{\sum \min(a_i, b_i)}{\sum a_i + \sum b_i - \sum \min(a_i, b_i)}}
#' 
#' \strong{Cosine Similarity:}
#' \deqn{cosine = \frac{\sum a_i b_i}{\sqrt{(\sum a_i^2)(\sum b_i^2)}}}
#' 
#' Where:
#' \itemize{
#'   \item{\eqn{a_i} and \eqn{b_i} are the abundances of species \eqn{i} in groups A and B, respectively.}
#' }
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' clonalOverlap(combined, 
#'               cloneCall = "aa", 
#'               method = "jaccard")
#'
#' @param input.data The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa),
#' VDJC gene + CDR3 nucleotide (strict) or a custom variable in the data. 
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the "overlap", "morisita", 
#' "jaccard", "cosine" indices or "raw" for the base numbers.
#' @param group.by The variable to use for grouping.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @param pcoa Plot PCoA plot instead of a heatmap.
#' @param pcoa.group.by grouping variable in PCoA plot, default to NULL.
#' @param point.size point size in PCoA plot.
#' @importFrom stringr str_sort str_to_title
#' @importFrom reshape2 melt
#' @importFrom stats quantile
#' @importFrom ape pcoa
#' @export
#' @concept Visualizing_Clones
#' @return ggplot of the overlap of clones by group
clonalOverlap <- function(input.data, 
                          cloneCall = "strict", 
                          method = NULL, 
                          chain = "both", 
                          group.by = NULL,
                          exportTable = FALSE,
                          palette = "inferno",
                          pcoa = FALSE,
                          pcoa.group.by = NULL,
                          point.size = 3){
    if(method == "morisita") {
      return_type = "freq"
    } else {
      return_type = "unique"
    }
  input.data <- .data.wrangle(input.data, 
                              group.by, 
                              .theCall(input.data, cloneCall, check.df = FALSE), 
                              chain)
  cloneCall <- .theCall(input.data, cloneCall)
    
    input.data <- input.data[order(names(input.data))]
    values <- str_sort(as.character(unique(names(input.data))), numeric = TRUE)
    input.data <- input.data[values]
    num_samples <- length(input.data[])
    names_samples <- names(input.data)
    length <- seq_len(num_samples)
    
    if (chain != "both") {
      for (i in seq_along(input.data)) {
        input.data[[i]] <- .off.the.chain(input.data[[i]], chain, cloneCall)
      }
    }
    
    #Selecting Index Function
    indexFunc <- switch(method,
                        "morisita" = .morisitaCalc,
                        "jaccard"  = .jaccardCalc,
                        "raw"      = .rawCalc,
                        "overlap"  = .overlapCalc,
                        "cosine"  = .cosineCalc,
                        stop("Invalid method provided"))
    
    #Calculating Index 
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    coef_matrix <- .calculateIndex(input.data, length, cloneCall, coef_matrix, indexFunc, return_type)
    
    #Data manipulation
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples

    if (exportTable == TRUE) { 
      return(coef_matrix) 
    }

    if (pcoa) {
        if (!is.null(pcoa.group.by)) {
            group <- unlist(lapply(input.data,
                function(x) unique(x[[pcoa.group.by]])))
        } else {
            group <- NULL
        }
        m <- as.matrix(coef_matrix)
        m[lower.tri(m)] <- t(m)[lower.tri(m)]
        m <- as.dist(1-m, upper=TRUE)
        res_pcoa <- ape::pcoa(m, correction="lingoes")
        plot <- .returnPcoa(res_pcoa, group=group, palette=palette,
            point.size=point.size, pcoa.group.by=pcoa.group.by)
        return(plot)
    }

    mat_melt <- suppressMessages(melt(as.matrix(coef_matrix)))
    
    mean_value <- mean(na.omit(mat_melt[,"value"]))
    
    plot <- ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
                geom_tile() + 
                geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], fill = NA, lwd= 0.25, color = "black") +
                labs(fill = str_to_title(method)) +
                geom_text(aes(label = round(value, digits = 3), 
                              color = ifelse(value <= mean_value,
                                             "white", "black")), 
                          na.rm = TRUE) +
                scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
                scale_color_identity() +
                theme_classic() + 
                theme(axis.title = element_blank())
    return(plot) 
}

# Helper function to prepare data
.prepareDataFrame <- function(df, cloneCall, return_type = "unique") {
  if (return_type == "unique") {
    return(unique(df[, cloneCall]))
  } else if (return_type == "freq") {
    temp_df <- data.frame(table(df[, cloneCall]))
    colnames(temp_df) <- c(cloneCall, 'Count')
    temp_df[, 2] <- as.numeric(temp_df[, 2])
    return(temp_df)
  }
}

# Helper function for common loop and conditional structure
.calculateIndex <- function(df, length, cloneCall, coef_matrix, indexFunc, return_type = "unique") {
  for (i in seq_along(length)) {
    df_i <- .prepareDataFrame(df[[i]], cloneCall, return_type)
    for (j in seq_along(length)) {
      if (i >= j) { next }
      df_j <- .prepareDataFrame(df[[j]], cloneCall, return_type)
      coef_matrix[i, j] <- indexFunc(df_i, df_j)
    }
  }
  return(coef_matrix)
}

# Morisita Index calculation function
.morisitaCalc <- function(df_i, df_j) {
  merged <- merge(df_i, df_j, by = names(df_i)[1], all = TRUE)
  merged[is.na(merged)] <- 0
  
  X <- sum(merged[, 2])
  Y <- sum(merged[, 3])
  
  num <- 2 * sum(merged[, 2] * merged[, 3])
  den <- ((sum(df_i[, 2]^2) / (X^2)) + (sum(df_j[, 2]^2) / (Y^2))) * X * Y
  
  return(num / den)
}

# Jaccard Index calculation function
.jaccardCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / (length(df_i) + length(df_j) - overlap))
}

# Raw Index calculation function
.rawCalc <- function(df_i, df_j) {
  return(length(intersect(df_i, df_j)))
}

# Overlap Index calculation function
.overlapCalc <- function(df_i, df_j) {
  overlap <- length(intersect(df_i, df_j))
  return(overlap / min(length(df_i), length(df_j)))
}

# Overlap Index calculation function
.cosineCalc <- function(df_i, df_j) {
  all_species <- unique(c(df_i, df_j))
  vector_location1 <- as.integer(all_species %in% df_i)
  vector_location2 <- as.integer(all_species %in% df_j)
  
  return(sum(vector_location1 * vector_location2) / 
           (sqrt(sum(vector_location1^2)) * sqrt(sum(vector_location2^2))))
}
