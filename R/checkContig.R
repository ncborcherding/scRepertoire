#' Assess the number of NAs in clonotypes
#'
#' This function takes the output from combineTCR(), combineBCR(), 
#' combineExpression() or expression2List() and quantifies the 
#' percent of NA values in either the sequence information of 
#' the VDJC gene usage. If a data frame output for the raw numbers 
#' and percentage of NAs are preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' checkContig(combined, chain = "TRA", examine = "chain", group.by = "sample")
#'
#' @param df The product of combineTCR(), combineBCR(), expression2List(), or combineExpression().
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param examine indicate either "chain" for sequence-based NA or 
#' "gene" for NA percentages across the VDJC genes.
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param exportTable Returns the data frame used for forming the graph
#' @import dplyr
#' @importFrom stringr str_split
#' @export
#' @return ggplot of the percentage of NA values
#' 
checkContig <- function(df, 
                         chain = "TRA", 
                         examine = "gene",
                         group.by = NULL, 
                         split.by = NULL, 
                         exportTable = FALSE) {
  
  x.axis <- c(split.by, group.by) #condensing for ploting.
  x.axis <- x.axis[1]
  if(is.null(x.axis)) {
    stop("Please select a group.by or split.by variable.")
  }
  if(chain %in% c("TRA", "TRG", "IGL")) {
    genes <- c("v_gene", "j_gene", "c_gene")
  } else {
    genes <- c("v_gene", "d_genes", "j_gene", "c_gene")
  }
  df <- list.input.return(df,split.by)
  df <- checkBlanks(df, "CTaa")
  for (i in seq_along(df)) {
    df[[i]] <- off.the.chain(df[[i]], chain, "CTaa")
    df[[i]] <- off.the.chain(df[[i]], chain, "CTnt")
    df[[i]] <- off.the.chain(df[[i]], chain, "CTgene")
    df[[i]][df[[i]] == "NA"] <- NA
  }
  
  names <- names(df)
  
  dat <- bind_rows(df)
  if (examine == "chain") {
    if(!is.null(split.by)) {
      na.by.chain <- dat %>%
        group_by(dat[,split.by]) %>%
        summarise(NA.values = sum(is.na(CTaa)), 
                  total.cells = n())
      colnames(na.by.chain)[1] <- split.by
    } else if(!is.null(group.by)) {
      na.by.chain <- dat %>%
        group_by(dat[,group.by]) %>%
        summarise(NA.values = sum(is.na(CTaa)), 
                  total.cells = n())
      colnames(na.by.chain)[1] <- group.by
    } else {
      stop("Please select a split.by or group.by variable")
    }
    if(exportTable) {
      return(na.by.chain)
    }
    na.by.chain <- na.by.chain %>%
      mutate(value = round((NA.values/total.cells)*100, 2)) %>%
      as.data.frame()
    
    na.by.chain$variable<- "chain"
    plot.data <- na.by.chain
  } else if (examine == "gene") {
    out <- str_split(dat[,"CTgene"], ";", simplify = TRUE)[,1]
    out <- data.frame(str_split(out, "[.]", simplify = TRUE))
    out [out ==""] <- NA
    out [out =="None"] <- NA
    dat <- cbind(dat[,c(x.axis)], out)
    colnames(dat)[1] <- x.axis
    #############
    #Summarise the NAs by column/group
    na.by.gene <- dat %>%
      group_by(dat[,x.axis]) %>%
      summarise(across(colnames(dat)[seq_len(ncol(dat))[-1]], list(~ sum(is.na(.)))))
    #Summarise number by column/group
    colnames(na.by.gene) <- c(x.axis, genes)
    total.value  <- dat %>%
      group_by(dat[,x.axis]) %>% 
      summarise(total.cells = n())
    na.by.gene <- data.frame(na.by.gene, total = total.value[,ncol(total.value)])
    if(exportTable) {
      return(na.by.gene)
    }
    na.by.gene <- na.by.gene %>%
      mutate(across(genes, ~round((./total.cells)*100,2)))
    na.by.gene <- melt(na.by.gene[,seq_len(length(genes) + 1)], 
                       id.vars = x.axis)
    plot.data <- na.by.gene
    
  }
  col <- length(unique(plot.data[,x.axis]))
  plot <- ggplot(plot.data, aes(x = variable, y = value)) + 
    geom_boxplot(width = 0.5, fill = "grey", outlier.alpha = 0) + 
    geom_jitter(aes(color = plot.data[,x.axis]), size = 3) + 
    labs(color = x.axis) + 
    theme_classic() + 
    ylab("Percent of NAs") + 
    scale_color_manual(values = colorblind_vector(col)) + 
    theme(axis.title.x = element_blank())
  return(plot)
}
