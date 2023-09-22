#' Examining the clonal overlap between groups or samples
#'
#' This functions allows for the calculation and visualizations of the 
#' overlap coefficient, morisita, or jaccard index for clonotypes 
#' using the product of combineTCR(), combineBCR() or expression2list(). 
#' The overlap coefficient is calculated using the intersection of clonotypes 
#' divided by the length of the smallest component. 
#' If a matrix output for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' clonalOverlap(combined, 
#'               cloneCall = "gene", 
#'               method = "jaccard")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the overlap, 
#'  "overlap" 
#   "morisita", "jaccard" indices, "cosine" similarity or "raw" 
#' for the base numbers.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @importFrom stringr str_sort str_to_title
#' @importFrom reshape2 melt
#' @importFrom stats quantile
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list
clonalOverlap <- function(df, 
                          cloneCall = "strict", 
                          method = NULL, 
                          chain = "both", 
                          split.by = NULL,
                          exportTable = FALSE,
                          palette = "inferno"){
    if(method == "morisita") {
      return_type = "freq"
    } else {
      return_type = "unique"
    }
    df <- list.input.return(df, split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    df <- df[order(names(df))]
    values <- str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]
    num_samples <- length(df[])
    names_samples <- names(df)
    length <- seq_len(num_samples)
    
    if (chain != "both") {
      for (i in seq_along(df)) {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
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
    coef_matrix <- .calculateIndex(df, length, cloneCall, coef_matrix, indexFunc, return_type)
    
    #Data manipulation
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples

    if (exportTable == TRUE) { 
      return(coef_matrix) 
    }
    mat_melt <- suppressMessages(melt(as.matrix(coef_matrix)))
    
    mean_value <- mean(na.omit(mat_melt[,"value"]))
    
    plot <- ggplot(mat_melt, aes(x=Var1, y=Var2, fill=value)) +
                geom_tile() + 
                geom_tile(data = mat_melt[!is.na(mat_melt[,"value"]),], fill = NA, lwd= 0.25, color = "black") +
                labs(fill = str_to_title(method)) +
                geom_text(aes(label = round(value, digits = 3), 
                              color = ifelse(value <= mean_value,
                                             "white", "black"))) +
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