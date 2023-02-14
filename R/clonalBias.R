#' Examine clonotype bias
#' 
#' Clonotype bias method was developed and outlined from a single-cell 
#' \href{https://www.biorxiv.org/content/10.1101/2021.09.20.458613v1}{manuscript} 
#' characterizing CD4 responses to acute and chronic infection. The metric seeks to 
#' quantify how individual clones are skewed towards a specific cellular 
#' compartment or cluster. A clonotype bias of 1 indicates that a clonotype 
#' is composed of cells from a single compartment or cluster, while a clonotype
#' bias of 0 matches the background subtype distribution. 
#' 
#' @examples
#' \dontrun{
#' Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' screp_example <- get(data("screp_example"))
#' sce <- suppressMessages(Seurat::UpdateSeuratObject(screp_example))
#' 
#' #Using combineExpresion()
#' sce <- combineExpression(combined, sce)
#' 
#' #Using occupiedscRepertoire()
#' clonotypeBias(sce, cloneCall = "CTaa", split.by = "Patient", group.by = "seurat_clusters",
#' n.boots = 20, min.expand = 2)
#' }
#' 
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param group.by The column header used for comparisons of bias.
#' @param split.by The column header used for calculating the baseline frequencies.
#' For example, "Type" for tumor vs peripheral blood comparison 
#' @param n.boots number of bootstraps to downsample
#' @param min.expand clonotype frequency cut off for the purpose of comparison
#' @param exportTable Returns the data frame used for forming the graph
#' @import ggplot2
#' @importFrom stringr str_sort
#' @export
#' @return Returns ggplot of the clonotype bias
clonotypeBias <- function(df, 
                          cloneCall="strict", 
                          split.by=NULL, 
                          group.by=NULL, 
                          n.boots = 20,
                          min.expand=10,
                          exportTable = FALSE) {
  
  bias <- get_clono_bias(df, 
                         split.by = split.by, 
                         group.by = group.by , 
                         cloneCall=cloneCall, 
                         min.expand=min.expand)
  df_shuffle.list <- list()
  for (ii in seq_len(n.boots)) {
    df_shuffle.list[[ii]] <- get_clono_bias(df, split.by = split.by,
                                group.by = group.by, 
                                cloneCall=cloneCall, 
                                min.expand=min.expand, 
                                do.shuffle = T, 
                                seed=ii)
  }
  df_shuffle <- Reduce(rbind, df_shuffle.list)
  
  stat.summary <- df_shuffle %>%
                      group_by(ncells) %>%
                      summarise(mean = mean(bias), 
                                std = sd(bias))
  
  corrected_p <- 1-(0.05/nrow(bias))
  bias$Top_state <- factor(bias$Top_state, str_sort(unique(bias$Top_state), numeric = TRUE))
  
  bias$Z.score <- NA
  for(i in seq_len(nrow(bias))) {
      stat.pos <- bias[i,]$ncells
      row.pull <- stat.summary[stat.summary$ncells == stat.pos,]
      z.score <- (bias[i,]$bias - row.pull$mean)/row.pull$std
      bias$Z.score[i] <- z.score
  }
  
  plot <- ggplot(bias, aes(x=ncells,y=bias)) + 
    geom_point(aes(colour=Top_state)) + 
    quiet(stat_quantile(data=df_shuffle, quantiles = c(corrected_p), method = "rqss", lambda=3)) + 
    theme_classic() + 
    xlab("Clone Size") + 
    ylab("Clonotype Bias")
  if (exportTable == TRUE) { 
    return(bias) 
  }
  return(plot) 
}

get_clono_bg <- function(df, 
                         split.by=split.by, 
                         group.by=group.by, 
                         cloneCall=cloneCall, 
                         min.expand=10) {
  if (!is.null(split.by)) {
    df <- list.input.return(df, split.by)
  } else {
    if (inherits(x=df, what ="Seurat") | inherits(x=df, what ="SummarizedExperiment")) {
      df <- list("Object" = grabMeta(df))
    } 
  }
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  
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

#Code Derived from 
get_clono_bias <- function(df, 
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
  
  bg <- get_clono_bg(df, split.by=split.by, 
                     min.expand = min.expand, 
                     group.by = group.by, 
                     cloneCall = cloneCall)
  if (!is.null(split.by)) {
    df <- list.input.return(df, split.by)
  } else {
    if (inherits(x=df, what ="Seurat") | inherits(x=df, what ="SummarizedExperiment")) {
      df <- list("Object" = grabMeta(df))
    } 
  }
  cloneCall <- theCall(cloneCall)
  df <- checkBlanks(df, cloneCall)
  
  for (s in names(bg)) {
    
    #All cells have the same type. Cannot calculate bias for this sample
    if (length(bg[[s]]) > 1) {  
      
      clones <-  table(df[[s]][,cloneCall])
      clones <- clones[clones>=min.expand]
      subtypes <- names(bg[[s]])
      
      if (length(clones)>0) {
        clones <- sort(clones, decreasing = T)
        expanded <- df[[s]][which(df[[s]][,cloneCall] %in% names(clones)), c(group.by, cloneCall)]
        
        if (do.shuffle) {  #reshuffle annotation column
          set.seed(seed)
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
          dat[nrow(dat) + 1,] = new_row
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
  