#' Examine clonotype bias
#' 
#' The metric seeks to quantify how individual clones are skewed towards 
#' a specific cellular compartment or cluster. A clonotype bias of *1* - 
#' indicates that a clonotype is composed of cells from a single 
#' compartment or cluster, while a clonotype bias of *0* - matches the 
#' background subtype distribution. Please read and cite the following
#' \href{https://pubmed.ncbi.nlm.nih.gov/35829695/}{manuscript} 
#' if using \code{\link{clonalBias}}.
#' 
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' #Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' scRep_example$Patient <- substring(scRep_example$orig.ident,1,3)
#' 
#' #Using clonalBias()
#' clonalBias(scRep_example, 
#'               cloneCall = "aa", 
#'               split.by = "Patient", 
#'               group.by = "seurat_clusters",
#'               n.boots = 5, 
#'               min.expand = 2)
#' 
#' 
#' @param df The product of \code{\link{combineExpression}}.
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param group.by The variable to use for calculating bias
#' @param split.by The variable to use for calculating the baseline frequencies.
#' For example, "Type" for lung vs peripheral blood comparison 
#' @param n.boots number of bootstraps to downsample.
#' @param min.expand clonotype frequency cut off for the purpose of comparison.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom stringr str_sort
#' @export
#' @return Returns ggplot of the clonotype bias
clonalBias <- function(df, 
                       cloneCall="strict", 
                       split.by=NULL, 
                       group.by=NULL, 
                       n.boots = 20,
                       min.expand=10,
                       exportTable = FALSE, 
                       palette = "inferno") {
  .checkSingleObject(df)
  cloneCall <- .theCall(cloneCall)
  #Calculating bias
  bias <- .get_clono_bias(df, 
                         split.by = split.by, 
                         group.by = group.by , 
                         cloneCall=cloneCall, 
                         min.expand=min.expand)
  df_shuffle.list <- list()
  #Bootstrapping
  for (ii in seq_len(n.boots)) {
    df_shuffle.list[[ii]] <- .get_clono_bias(df, split.by = split.by,
                                            group.by = group.by, 
                                            cloneCall=cloneCall, 
                                            min.expand=min.expand, 
                                            do.shuffle = TRUE, 
                                            seed=ii)
  }
  df_shuffle <- Reduce(rbind, df_shuffle.list)
  
  #Sumarrising boot straps
  stat.summary <- df_shuffle %>%
    group_by(ncells) %>%
    summarise(mean = mean(bias), 
              std = sd(bias))
  
  corrected_p <- 1-(0.05/nrow(bias))
  bias$Top_state <- factor(bias$Top_state, str_sort(unique(bias$Top_state), numeric = TRUE))
  
  #Calculating Bias Z-score
  bias$Z.score <- NA
  for(i in seq_len(nrow(bias))) {
    stat.pos <- bias[i,]$ncells
    row.pull <- stat.summary[stat.summary$ncells == stat.pos,]
    z.score <- (bias[i,]$bias - row.pull$mean)/row.pull$std
    bias$Z.score[i] <- z.score
  }
  
  #Attaching the cloneSize of original combineExpression()
  meta <- .grabMeta(df)
  meta <- meta[meta[,cloneCall] %in% bias[,"Clone"],]
  meta <- unique(meta[,c(cloneCall, split.by, "cloneSize")])
  
  bias$cloneSize <- NA
  for(i in seq_len(nrow(meta))) {
    split <- bias[,1][i]
    clone <- bias[,3][i]
    bias$cloneSize[i] <- as.vector(meta[which(meta[,cloneCall] == clone & meta[,split.by] == split),"cloneSize"])
  }
  
  bias$cloneSize <- factor(bias$cloneSize , rev(levels(meta[,"cloneSize"])))
  #Plotting 
  plot <- ggplot(bias, aes(x=ncells,y=bias)) + 
    geom_point(aes(fill=Top_state, size = cloneSize), shape = 21, stroke = 0.25) + 
    .quiet(stat_quantile(data=df_shuffle, 
                         quantiles = c(corrected_p), 
                         method = "rqss", 
                         lambda=3, 
                         color = "black", 
                         lty = 2)) + 
    scale_fill_manual(values = .colorizer(palette,  length(unique(bias[,"Top_state"])))) + 
    theme_classic() + 
    xlab("Clone Size") + 
    ylab("Clonotype Bias")
  if (exportTable == TRUE) { 
    return(bias) 
  }
  return(plot) 
}

#Background summary of clones
.get_clono_bg <- function(df, 
                         split.by=split.by, 
                         group.by=group.by, 
                         cloneCall=cloneCall, 
                         min.expand=10) {
  df <- .data.wrangle(df, split.by, cloneCall, "both")
  
  bg <- list()
  for (s in seq_along(df)) {
    
    clones <-  table(df[[s]][,cloneCall])
    clones <- clones[clones>=min.expand]
    
    if (length(clones)>0) {
      
      expanded <- df[[s]][which(df[[s]][,cloneCall] %in% names(clones)), group.by]
      bg[[s]] <- table(expanded) / sum(table(expanded))
    }
  }
  names(bg) <- names(df)
  return(bg)
}

#Clone Bias
.get_clono_bias <- function(df, 
                           split.by=NULL, 
                           group.by=NULL, 
                           cloneCall=cloneCall, 
                           min.expand=10,
                           do.shuffle=FALSE, 
                           seed=123) {
  
  dat <- data.frame(Sample=character(),
                    Clone_i=character(),
                    Clone=character(),
                    ncells=integer(),
                    Top_state=double(), 
                    freq=double(),
                    freq_diff=double(),
                    bias=double()) 
  
  bg <- .get_clono_bg(df, 
                     split.by=split.by, 
                     min.expand = min.expand, 
                     group.by = group.by, 
                     cloneCall = cloneCall)
  if (!is.null(split.by)) {
    df <- .list.input.return(df, split.by)
  } else {
    if (inherits(x=df, what ="Seurat") | inherits(x=df, what ="SummarizedExperiment")) {
      df <- list("Object" = .grabMeta(df))
    } 
  }
  cloneCall <- .theCall(cloneCall)
  df <- .checkBlanks(df, cloneCall)
  
  for (s in names(bg)) {
    
    #All cells have the same type. Cannot calculate bias for this sample
    if (length(bg[[s]]) > 1) {  
      
      clones <-  table(df[[s]][,cloneCall])
      clones <- clones[clones>=min.expand]
      subtypes <- names(bg[[s]])
      
      if (length(clones)>0) {
        clones <- sort(clones, decreasing = TRUE)
        expanded <- df[[s]][which(df[[s]][,cloneCall] %in% names(clones)), c(group.by, cloneCall)]
        
        if (do.shuffle) {  #reshuffle annotation column
          expanded[[group.by]] <- sample(expanded[[group.by]])
        }
        
        for (i in seq_along(clones)) {
          this <- names(clones)[i]
          this.s <- paste0(this, "_", s)
          
          ncells <- clones[i]
          sub <- expanded[which(expanded[, cloneCall] == this), group.by]   
          sub <- factor(sub, levels=subtypes)
          
          comp <- table(sub) / sum(table(sub))
          diff <- comp - bg[[s]]
          diff_norm <- round((comp - bg[[s]])/(1-bg[[s]]), 3)
          top_state <- names(which.max(diff_norm))[1]
          
          new_row <- c(s, this.s, this, as.integer(ncells), top_state, 
                       as.double(comp[top_state]), 
                       as.double(diff[top_state]), 
                       as.double(diff_norm[top_state]))
          dat[nrow(dat) + 1,] <- new_row
        }
      }
    }
  }
  dat$ncells <- as.numeric(dat$ncells)
  dat$freq <- as.double(dat$freq)
  dat$freq_diff <- as.double(dat$freq_diff)
  dat$bias <- as.double(dat$bias)
  return(dat)
}
