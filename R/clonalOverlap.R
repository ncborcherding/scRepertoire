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
#' clonalOverlap(combined, cloneCall = "gene", method = "overlap")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The method to calculate the overlap, either the "overlap" 
#' coefficient, "morisita", "jaccard" indices, "cosine" similarity or "raw" 
#' for the base numbers.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @export
#' @return ggplot of the clonotypic overlap between elements of a list
clonalOverlap <- function(df, 
                          cloneCall = "strict", 
                          method = NULL, 
                          chain = "both", 
                          split.by = NULL, 
                          exportTable = FALSE,
                          palette = "inferno"){
    df <- list.input.return(df, split.by)
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    df <- df[order(names(df))]
    values <- str_sort(as.character(unique(names(df))), numeric = TRUE)
    df <- df[quiet(dput(values))]
    num_samples <- length(df[])
    names_samples <- names(df)
    coef_matrix <- data.frame(matrix(NA, num_samples, num_samples))
    colnames(coef_matrix) <- names_samples
    rownames(coef_matrix) <- names_samples
    length <- seq_len(num_samples)
    if (chain != "both") {
      for (i in seq_along(df)) {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
    }
    
    coef_matrix <- switch(method,
                     "cosine" = .cosineIndex(df, length, cloneCall, coef_matrix),
                     "jaccard" = .jaccardIndex(df, length, cloneCall, coef_matrix),
                     "morisita" = .morisitaIndex(df, length, cloneCall, coef_matrix),
                     "overlap" = .overlapIndex(df, length, cloneCall, coef_matrix),
                     "raw" = .rawIndex(df, length, cloneCall, coef_matrix),
                     "Invalid method specified")
    coef_matrix <- as.data.frame(coef_matrix)
    coef_matrix$names <- rownames(coef_matrix)
    if (exportTable == TRUE) { 
      return(coef_matrix) 
    }
    coef_matrix <- suppressMessages(melt(coef_matrix))
    coef_matrix$variable <- factor(coef_matrix$variable, levels = values)
    coef_matrix$names <- factor(coef_matrix$names, levels = values)
    
    tertile_values <- quantile(na.omit(coef_matrix[,"value"]), probs = c(1/3,2/3))
    
    plot <- ggplot(coef_matrix, aes(x=names, y=variable, fill=value)) +
                geom_tile() + 
                geom_tile(data = coef_matrix[!is.na(coef_matrix[,"value"]),], fill = NA, lwd= 0.25, color = "black") +
                labs(fill = method) +
                geom_text(aes(label = round(value, digits = 3), 
                              color = ifelse(value <= as.vector(tertile_values[1]),
                                             "white", "black"))) +
                scale_fill_gradientn(colors = .colorizer(palette, 7), na.value = "white") +
                scale_color_identity() +
                theme_classic() + 
                theme(axis.title = element_blank())
    return(plot) 
}

#Calculate the Morisita Index for Overlap Analysis
#' @author Massimo Andreatta, Nick Borcherding
.morisitaIndex <- function(df, length, cloneCall, coef_matrix) {
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- data.frame(table(df.i[,cloneCall]))
    colnames(df.i) <- c(cloneCall, 'Count')
    df.i[,2] <- as.numeric(df.i[,2])
    for (j in seq_along(length)){
      if (i >= j){ next }
      else { df.j <- df[[j]]
      df.j <- data.frame(table(df.j[,cloneCall]))
      colnames(df.j) <- c(cloneCall, 'Count')
      df.j[,2] <- as.numeric(df.j[,2])
      merged <- merge(df.i, df.j, by = cloneCall, all = TRUE)
      merged[is.na(merged)] <- 0
      X <- sum(merged[,2])
      Y <- sum(merged[,3])
      sum.df.i <- sum(df.i[,2]^2)
      sum.df.j <- sum(df.j[,2]^2)
      
      num <- 2 * sum(merged[, 2] * merged[, 3])
      den <- ((sum.df.i / (X^2) + sum.df.j / (Y^2)) * X * Y)
      
      coef.i.j <- num/den
      coef_matrix[i,j] <- coef.i.j
      }
    }
  }
  return(coef_matrix)
}


#Calculate the Jaccard Similarity Index for Overlap Analysis
.jaccardIndex <- function(df, length, cloneCall, coef_matrix) {
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,cloneCall]
    df.i_unique <- unique(df.i)
    for (j in seq_along(length)){
      if (i >= j){ next }
      
      df.j <- df[[j]]
      df.j <- df.j[,cloneCall]
      df.j_unique <- unique(df.j)
      overlap <- length(intersect(df.i_unique, 
                                  df.j_unique))
      coef_matrix[i,j] <- 
        overlap/(sum(length(df.i_unique), 
                     length(df.j_unique))-overlap)
    }
  }
  return(coef_matrix)
}

.rawIndex <- function(df, length, cloneCall, coef_matrix) {
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,cloneCall]
    df.i_unique <- unique(df.i)
    for (j in seq_along(length)){
      if (i >= j){ next }
      df.j <- df[[j]]
      df.j <- df.j[,cloneCall]
      df.j_unique <- unique(df.j)
      overlap <- length(intersect(df.i_unique, 
                                  df.j_unique))
      coef_matrix[i,j] <- overlap
    }
  }
  return(coef_matrix)
}


#Calculate the Overlap Coefficient for Overlap Analysis
#' @author Nick Bormann, Nick Borcherding
.overlapIndex <- function(df, length, cloneCall, coef_matrix) {
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,c(cloneCall)]
    df.i_unique <- unique(df.i)
    for (j in seq_along(length)){
      if (i >= j){ next }
      else { df.j <- df[[j]]
      df.j <- df.j[,c(cloneCall)]
      df.j_unique <- unique(df.j)
      overlap <- length(intersect(df.i_unique, 
                                  df.j_unique))
      coef_matrix[i,j] <- 
        overlap/min(length(df.i_unique), 
                    length(df.j_unique)) } } }
  return(coef_matrix)
}

.cosineIndex <- function(df, length, cloneCall, coef_matrix) {
  for (i in seq_along(length)){
    df.i <- df[[i]]
    df.i <- df.i[,cloneCall]
    df.i_unique <- unique(df.i)
    for (j in seq_along(length)){
      if (i >= j){ next }
      else { 
        df.j <- df[[j]]
        df.j <- df.j[,cloneCall]
        df.j_unique <- unique(df.j)
        all_species <- unique(c(df.j_unique, df.j_unique))
        vector_location1 <- as.integer(all_species %in% df.i_unique)
        vector_location2 <- as.integer(all_species %in% df.j_unique)
        
        coef_matrix[i,j] <- 
          sum(vector_location1 * vector_location2) / (sqrt(sum(vector_location1^2)) * sqrt(sum(vector_location2^2)))
      } 
    } 
  }
  coef_matrix <- as.matrix(coef_matrix)
  coef_matrix[is.nan(coef_matrix)] <- 0
  return(coef_matrix)
}

