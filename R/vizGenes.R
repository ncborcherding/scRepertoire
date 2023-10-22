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
#' vizGenes(combined, 
#'          x.axis = "TRBV",
#'          y.axis = NULL,
#'          plot = "heatmap")
#'
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param plot The type of plot to return - heatmap or barplot 
#' @param x.axis Gene segments to separate the x-axis, such as "TRAV", "TRBD", "IGKJ".
#' @param y.axis Variable to separate the y-axis, can be both categorical or other gene 
#' gene segments, such as "TRAV", "TRBD", "IGKJ".
#' @param group.by Variable in which to group the diversity calculation.
#' @param order Categorical variable to organize the x-axis, either "gene" or "variance"
#' @param scale Converts the individual count of genes to proportion using the total 
#' respective repertoire size 
#' @param exportTable Returns the data frame used for forming the graph.
#' @param palette Colors to use in visualization - input any \link[grDevices]{hcl.pals}.
#' @import ggplot2
#' @importFrom stringr str_split
#' @importFrom stats sd
#' @importFrom dplyr bind_rows
#' @export
#' @concept Visualizing_Clones
#' @return ggplot bar diagram or heatmap of gene usage

vizGenes <- function(df, 
                     x.axis = "TRBV",
                     y.axis = NULL,
                     group.by = NULL, 
                     plot = "heatmap",  
                     order = "gene",
                     scale = TRUE, 
                     exportTable = FALSE,
                     palette = "inferno") {
  element.names <- NULL
  sco <- is_seurat_object(df) | is_se_object(df)
  
  #Extracting group.by in case of null
  if(!grepl("TRA|TRB|TRG|TRD|IGH|IGL|IGK", y.axis) && !is.null(y.axis)) {
    group.by <- y.axis
  }
  #If single-cell object, need group.by
  if(is.null(group.by) & sco) {
    group.by <- "ident"
  }
  
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
  if (any(y.axis %!in% colnames(df)) | !is.null(y.axis)) {
    if (grepl("TRA|TRB|TRG|TRD|IGH|IGL|IGK", y.axis)) {
      df <- .select.gene(df, y.axis, y.axis)
      colnames(df)[ncol(df)] <- y.axis ######## Check if need this
    } else {
      y.axis <- "element.names"
    }
  } else {
    y.axis <- "element.names"
  }
  
  #Filtering steps
  df[df == "NA"] <- NA
  df[df == ""] <- NA
  if(!is.null(y.axis) & any(is.na(df[,c(y.axis)]))) {
    df <- subset(df, !is.na(df[,c(y.axis)]))
  }
  if(!is.null(x.axis) & any(is.na(df[,c(x.axis)]))) {
    df <- subset(df, !is.na(df[,c(x.axis)]))
  }
  
  if (!scale) {
    col.lab <- "Total n"
    values <- "count"
  } else {
    col.lab <- "Scaled Values"
    values <- "proportion"
  }
  
  #Calculating the summary values
  if (!is.null(y.axis) && y.axis != "element.names") {
    mat <- df %>%
      group_by(element.names, df[,x.axis], df[,y.axis]) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      group_by(element.names) %>%
      mutate(total = sum(count), 
             proportion = count / total) %>%
      as.data.frame() %>%
      na.omit()
      colnames(mat)[2:3] <- c("x.axis", "y.axis")
      mat$n <- mat[,values]
      mat <- mat %>%
        group_by(y.axis, x.axis) %>%
        mutate(varcount = sum(n), 
               sd = sd(n, na.rm = TRUE),
               mean = mean(n)) %>%
        as.data.frame() 
    
  } else {
    mat <- df %>%
      group_by(element.names, df[,x.axis]) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      group_by(element.names) %>%
      mutate(total = sum(count), 
             proportion = count / total) %>%
      as.data.frame() %>%
      na.omit()
      colnames(mat)[1:2] <- c("y.axis", "x.axis")
      mat$n <- mat[,values]
      mat <- mat %>%
        group_by(y.axis, x.axis) %>%
        mutate(varcount = sum(n), 
               sd = sd(n, na.rm = TRUE),
               mean = mean(n)) %>%
        as.data.frame() 
  }
  
  
  if (order == "variance") {
    varOrder <- unique(mat[order(mat$varcount, decreasing = TRUE),"x.axis"])
  } else {
    varOrder <- str_sort(unique(mat$x.axis), numeric = TRUE)
  }
    mat[,"x.axis"] <- factor(mat[,"x.axis"], levels = varOrder)
  if (plot == "barplot") {
    mat2 <- unique(mat[,c("x.axis", "y.axis", "sd", "mean")])
    plot <- ggplot(mat2, aes(x=x.axis, y = mean)) + 
      geom_bar(stat = "identity", color = "black", lwd = 0.25) + 
      geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2,
                    position=position_dodge(.9)) + 
      theme_classic() + 
      theme(axis.title.x = element_blank(), #remove titles
            axis.title.y = element_blank(), 
            axis.ticks.x = element_blank(), #removes ticks
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5))) + 
      facet_grid(y.axis~.)
  } else if (plot == "heatmap") {
    
    plot <- ggplot(mat, aes(x=x.axis, y = y.axis)) + 
      geom_tile(aes(fill =log(n)), color = "black", lwd = 0.25) + 
      theme_classic() + 
      labs(fill = paste0("Log(",col.lab, ")")) + 
      theme(axis.title.x = element_blank(), #remove titles
            axis.ticks.x = element_blank(), #removes ticks
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5)), 
            axis.title.y = element_blank(), 
            axis.text.y = element_text(size=rel(0.5))) + 
      scale_fill_gradientn(colors = .colorizer(palette,5))
  }
  if (exportTable == TRUE) { 
    return(mat) 
  }
  return(plot)
  
  #Need to test:
  #Different chains - both heatmap and barplot
  #grouped elements - both heatmap and barplot
  
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
  df[,label] <- NA
  df[,label] <- genes
  return(df)
}
