#' Visualizing the distribution of gene usage
#' 
#' This function will allow for the visualizing the distribution 
#' of the any VDJ and C gene of the TCR or BCR using heatmap or
#' bar chart. This function requires assumes two chains were used in 
#' defining clonotype, if not, it will default to the only chain 
#' present regardless of the chain parameter.
#'
#' @examples
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' vizGenes(combined, gene = "V", chain = "TRB", plot = "bar", scale = TRUE)
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param gene Which part of the immune receptor to visualize - V, D, J, C
#' @param chain indicate the specific chain should be used - 
#' e.g. "TRA", "TRG", "IGH", "IGL" (no both option here)
#' @param plot The type of plot to return - heatmap or bar 
#' @param y.axis Variable to separate the y-axis, can be both categorical or other gene 
#' gene segments such as V, D, J, or C.
#' @param order Categorical variable to organize the x-axis, either "gene" or "variance"
#' @param scale Converts the individual count of genes to proportion using the total 
#' respective repertoire size 
#' @param group.by The variable to use for grouping.
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom stats sd
#' @importFrom dplyr bind_rows
#' @export
#' @return ggplot bar diagram or heatmap of gene usage

vizGenes <- function(df, 
                     gene = "V",
                     chain = "TRA", 
                     plot = "heatmap",  
                     y.axis = "sample", 
                     order = "gene",
                     scale = TRUE, 
                     group.by = NULL,
                     exportTable = FALSE,
                     palette = "inferno") {
  element.names <- NULL
  sco <- is_seurat_object(df) | is_se_object(df)
  df <- .list.input.return(df, split.by = group.by)
  if(!is.null(group.by) & !sco) {
    df <- .groupList(df, group.by)
  }
  for(i in seq_along(df)) {
    df[[i]] <- .off.the.chain(df[[i]], chain, "CTgene")
  }
  df <- .bound.input.return(df)
  if (y.axis %!in% colnames(df) | is.null(y.axis)) {
    if (y.axis %!in% c("V", "D", "J", "C")) {
      y.axis <- "element.names"
    } else {
      df <- .select.gene(df, chain, y.axis)
      colnames(df)[ncol(df)] <- y.axis
    }
  }
  df <- .select.gene(df, chain, gene)
  df <- subset(df, !is.na(df[,ncol(df)])) #remove NA values
  df <- subset(df, df[,ncol(df)] != "NA") #remove values that are character "NA"
  df <- subset(df, df[,ncol(df)] != "") #remove rows with non genes
  #df <- table(df[,ncol(df)], df[,y.axis])
  
  if (!is.null(y.axis) && y.axis != "element.names") {
    df <- df %>%
      group_by(df[,ncol(df)], df[,y.axis], element.names) %>%
      count() 
  } else {
    df <- df %>%
      group_by(df[,ncol(df)], element.names) %>%
      count() 
  }
  df <- df %>%
    group_by(element.names) %>%
    mutate(sum = sum(n))
  col.lab <- "Total n"
  if (scale == TRUE) {
    df[,"n"] <- df[,"n"]/df[,"sum"]
    col.lab <- "Scaled Values"
  } 
  colnames(df)[1:2] <- c("Var1", "Var2")
  df <- df %>%
    group_by(Var1, Var2) %>%
    mutate(varcount = sum(n), 
           sd = sd(n, na.rm = TRUE),
           mean = mean(n))
  if (order == "variance") {
    varOrder <- order(df$varcount, decreasing = TRUE)
    df$Var1 <- factor(df$Var1, levels = unique(df$Var1[varOrder]))
  }
  if (plot == "bar") {
    df2 <- unique(df[,c("Var1", "Var2", "sd", "mean")])
    plot <- ggplot(df2, aes(x=Var1, y = mean)) + 
      geom_bar(stat = "identity") + 
      geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                    position=position_dodge(.9)) + 
      theme_classic() + 
      theme(axis.title.x = element_blank(), #remove titles
            axis.title.y = element_blank(), 
            axis.ticks.x = element_blank(), #removes ticks
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5))) + 
      facet_grid(Var2~.)
  } else if (plot == "heatmap") {
    
    plot <- ggplot(df, aes(x=Var1, y = Var2)) + 
      geom_tile(aes(fill = n), color = "black") + 
      theme_classic() + 
      labs(fill = col.lab) + 
      theme(axis.title.x = element_blank(), #remove titles
            axis.ticks.x = element_blank(), #removes ticks
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5)), 
            axis.title.y = element_blank(), 
            axis.text.y = element_text(size=rel(0.5))) + 
      scale_fill_gradientn(colors = rev(.colorizer(palette,5)))
  }
  if (exportTable == TRUE) { return(df) }
  return(plot)
}


.vizGenes <- function(df, 
                     x.axis = "TRBV",
                     y.axis = "TRBJ",
                     plot = "heatmap",  
                     order = "gene",
                     scale = TRUE, 
                     group.by = NULL,
                     exportTable = FALSE,
                     palette = "inferno") {
  
  sco <- is_seurat_object(df) | is_se_object(df)
  df <- .list.input.return(df, split.by = group.by)
  if(!is.null(group.by) & !sco) {
    df <- .groupList(df, group.by)
  }
  
  df <- .bound.input.return(df)
  #Parsing x.axis if gene used
  if (x.axis %!in% colnames(df)) {
    if (grepl("TRA|TRB|TRG|TRD|IGH|IGL|IGK", x.axis)) {
      df <- .select.gene(df, x.axis, x.axis)
      colnames(df)[ncol(df)] <- x.axis ######## Check if need this
    }
  }
  
  #Parsing y.axis if gene used
  if (y.axis %!in% colnames(df) | is.null(y.axis)) {
    if (grepl("TRA|TRB|TRG|TRD|IGH|IGL|IGK", y.axis)) {
      df <- .select.gene(df, y.axis, y.axis)
      colnames(df)[ncol(df)] <- y.axis ######## Check if need this
    }
  }
  
  #Filtering steps
  df <- subset(df, !is.na(df[,c("x.axis", "y.axis")])) #remove NA values
  df <- subset(df, df[,ncol(df)] != "NA") #remove values that are character "NA"
  df <- subset(df, df[,ncol(df)] != "")
  
  
  ##################
  #Need to add the summary functions
  #Need to add the visualization functions
  
}

#Parsing the CTgene
.select.gene <- function(df, position, label) {
  chains <- c("TRAV" = 1, "TRBV" = 1, "TRGV" = 1, "TRDV" = 1, "IGHV" = 1, "IGLV" = 1, "IGKV" = 1,
              "TRAJ" = 2, "TRGJ" = 2, "IGKJ" = 2, "IKLJ" = 2,
              "TRBD" = 2, "TRDD" = 2, "IGHV" = 2, 
              "TRBJ" = 3, "TRDJ" = 2, "IGHJ" = 3)
  chain.poisiton <- chains[position]
  if(substring(position,3,3) %in% c("A", "G", "K", "L")) {
    chain.gene <- str_split(df[,"CTgene"], "_", simplify = TRUE)[,1]
  } else {
    chain.gene <- str_split(df[,"CTgene"], "_", simplify = TRUE)[,2]
  }
  genes <- str_split(chain.gene, "[.]", simplify = TRUE)[,chain.poisiton] 
  df[,label] <- genes
  return(df)
}