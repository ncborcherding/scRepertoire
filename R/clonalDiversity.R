#' Examine the clonal diversity of samples
#'
#' This function calculates traditional measures of diversity - Shannon, 
#' inverse Simpson, normalized entropy, Gini-Simpson, Chao1 index, and
#' abundance-based coverage estimators (ACE) measure of species evenness by sample or group. 
#' The function automatically down samples the diversity metrics using 
#' 100 boot straps The group parameter can be used to condense the individual 
#' samples. If a matrix output for the data is preferred, set exportTable = TRUE.
#' 
#' @details
#' The formulas for the indices and estimators are as follows:
#' 
#' \deqn{Shannon Index: H = - \sum p_i \cdot \log(p_i)}
#' 
#' \deqn{Inverse Simpson Index: D^{-1} = 1 / \sum p_i^2}
#' 
#' \deqn{Normalized Entropy: E^H = H / \log(S)}
#' 
#' \deqn{Gini-Simpson Index: 1 - D = 1 - \sum p_i^2}
#' 
#' \deqn{Chao1 Index: \hat{S}_{Chao1} = S + \left( \frac{n_1(n_1 - 1)}{2(n_2 + 1)} \right)}
#' 
#' \deqn{Abundance-based Coverage Estimator (ACE): \hat{S}_{ACE} = S_{abundant} + \frac{S_{rare}}{C_{rare}} + \left( \frac{S_{rare} - 1}{C_{rare}} \right) \cdot F_1}
#' 
#' Where:
#' \itemize{
#'   \item{ \( p_i \) is the proportion of species \( i \) in the dataset.}
#'   \item{ \( S \) is the total number of species.}
#'   \item{ \( n_1 \) and \( n_2 \) are the number of singletons and doubletons, respectively.}
#'   \item{ \( S_{abundant} \), \( S_{rare} \), \( C_{rare} \), and \( F_1 \) are parameters derived from the data.}
#' }
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' clonalDiversity(combined, cloneCall = "gene")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group.by Variable in which to group the diversity calculation
#' @param x.axis Additional variable in which to split the x.axis
#' @param group.by The variable to use for goruping
#' @param metrics The indices to use in diversity calculations - "shannon", "inv.simpson", 
#' "norm.entropy", "gini.simpson", "chao1", "ACE"
#' @param exportTable Exports a table of the data into the global environment 
#' in addition to the visualization
#' @param palette Colors to use in visualization - input any hcl.pals()
#' @param n.boots number of bootstraps to downsample in order to get mean diversity
#' @param return.boots export boot strapped values calculated - 
#' will automatically exportTable = TRUE
#' @param skip.boots remove downsampling and boot strapping from the calculation
#' @importFrom stringr str_sort str_split
#' @importFrom reshape2 melt
#' @importFrom dplyr sample_n
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
#' @author Andrew Malone, Nick Borcherding
clonalDiversity <- function(df, 
                            cloneCall = "strict", 
                            chain = "both",
                            group.by = NULL, 
                            x.axis = NULL, 
                            metrics = c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE"),
                            exportTable = FALSE, 
                            palette = "inferno",
                            n.boots = 100, 
                            return.boots = FALSE, 
                            skip.boots = FALSE) {
  if(return.boots) {
    exportTable <- TRUE
  }
  cloneCall <- .theCall(cloneCall)
  df <- .data.wrangle(df, group.by, cloneCall, chain)
  mat <- NULL
  sample <- c()
  if (!is.null(group.by) || !is.null(x.axis)) {
    df <- bind_rows(df, .id = "element.names")
    df$group.element <- paste0(df[,group.by], ".", df[,x.axis])
    #group.element.uniq <- unique(df$group.element)
    df <- split(df, f = df[,"group.element"])
  }
  min <- .short.check(df, cloneCall)
  for (i in seq_along(df)) {
      data <- as.data.frame(table(df[[i]][,cloneCall]))
      mat_a <- NULL
      sample <- c()
      if(skip.boots == TRUE) {
        sample <- .diversityCall(data)
        mat_a <- rbind(mat_a, sample)
        mat_a[is.na(mat_a)] <- 0
        mat <- rbind(mat, mat_a)
        mat <- as.data.frame(mat)
      } else {
        for (j in seq(seq_len(n.boots))) {
          x <- sample_n(data, min)
          sample <- .diversityCall(x)
          mat_a <- rbind(mat_a, sample)
        }
        mat_a[is.na(mat_a)] <- 0
        if(return.boots) {
          mat_a <- as.data.frame(mat_a)
          mat_a$sample <- names(df)[i]
          mat <- rbind(mat, mat_a)
        } else {
          mat_b<- colMeans(mat_a)
          mat_b<-as.data.frame(t(mat_b))
          mat <- rbind(mat, mat_b)
        }
      }
    }
    colnames(mat) <- c("shannon", "inv.simpson", "norm.entropy", "gini.simpson", "chao1", "ACE")
    mat <- mat[,colnames(mat) %in% metrics]
    if (!is.null(group.by)) {
      mat[,group.by] <- str_split(names(df), "[.]", simplify = TRUE)[,1]
    } else {
      group.by <- "Group"
      mat[,group.by] <- names(df)
    }
    if (!is.null(x.axis)) {
      mat[,x.axis] <- str_split(names(df), "[.]", simplify = TRUE)[,2]
    } else {
      x.axis <- "x.axis"
      mat[,x.axis] <- 1
    }
    if (exportTable == TRUE) { 
      return(mat) 
    }
    rownames(mat) <- names(df)
  
    mat_melt <- suppressMessages(melt(mat, id.vars = c(group.by, x.axis)))
    values <- str_sort(as.character(unique(mat_melt[,group.by])), 
                       numeric = TRUE)
    values <- .quiet(dput(values))
    mat_melt[,group.by] <- factor(mat_melt[,group.by], levels = values)
    if (x.axis == "x.axis") {
        plot <- ggplot(mat_melt, aes(x=1, y=as.numeric(value)))
    } else {
      plot <- ggplot(mat_melt, aes(x=mat_melt[,x.axis], y=as.numeric(value)))
    }
    plot <- plot +
      geom_boxplot(outlier.alpha = 0) +
      geom_jitter(aes(fill = mat_melt[,group.by]), size = 3, shape = 21, stroke = 0.25, color = "black") + 
      labs(fill = "Group") +
      ylab("Index Score") +
      scale_fill_manual(values = .colorizer(palette,length(unique(mat_melt[,group.by])))) +
    facet_wrap(~variable, scales = "free", ncol = length(metrics)) +
      theme_classic() + 
      theme(axis.title.x = element_blank())
    if (x.axis == "x.axis") { 
      plot <- plot + theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
      }
  
  return(plot) 
}


.diversityCall <- function(data) {
  shannon <- .shannon(data[,"Freq"])
  inv_simpson <- .invsimpson(data[,"Freq"])
  norm_entropy <- .normentropy(data[,"Freq"]) 
  gini_simpson <- .ginisimpson(data[,"Freq"]) 
  chao1 <- .chao1(data[,"Freq"])
  ACE <- .ACE(data[,"Freq"])
  out <- c(shannon, inv_simpson, norm_entropy, gini_simpson, chao1,ACE)
  return(out)
}


.shannon <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)))
}
.normentropy <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(-sum(p * log(p)) / log(length(p)))
}
.invsimpson <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 / sum(p^2))
}
.ginisimpson <- function(p){
  p <- p[which(p > 0)]
  p <- p / sum(p)
  p <- p[which(p > 0)]
  return(1 - sum(p^2))
}

.chao1 <- function(p){
  n1 <- sum(p == 1)
  n2 <- sum(p == 2)
  S_obs <- length(p)
  # Chao1 index calculation
  if(n1 > 1 && n2 > 0) {
    chao1 <- S_obs + (n1 * (n1 - 1)) / (2 * (n2 + 1))
  } else {
    # In cases where n1 <= 1 or n2 == 0, Chao1 is undefined
    chao1 <- NA
  }
  return(chao1)
}

.ACE <- function(p) {
  q = 10
  S_abund <- sum(p > q)
  rare_data <- p[p <= q]
  S_rare <- length(rare_data)
  n_rare <- sum(rare_data)
  
  # Calculate C_ACE
  C_ACE <- sum(p) / n_rare
  
  # Calculate gamma
  gamma <- 0
  for(i in seq_len(q)) {
    f_i <- sum(rare_data == i)
    gamma <- gamma + (1 - i / q)^f_i
  }
  
  # Calculate ACE
  ACE <- S_abund + (S_rare / C_ACE) + (1 - C_ACE) * gamma
  return(ACE)
}


