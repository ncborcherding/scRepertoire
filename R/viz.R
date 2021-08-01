#' Quantify the unique clonotypes in the filtered contigs.
#'
#' This function takes the output from combineTCR(), combineBCR(), or 
#' expression2List() and quantifies unique clonotypes. The unique clonotypes 
#' can be either reported as a raw output or scaled to the total number of 
#' clonotypes recovered using the scale parameter. Multiple sequencing 
#' runs can be group together using the group parameter. If a matrix output 
#' for the data is preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' quantContig(combined, cloneCall="gene+nt", scale = TRUE)
#'
#' @param df The product of combineTCR() combineBCR() or expression2List().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group The column header used for grouping.
#' @param scale Converts the graphs into percentage of unique clonotypes.
#' @param exportTable Returns the data frame used for forming the graph
#' @import ggplot2
#' @export
#' @return ggplot of the total or relative unique clonotypes
quantContig <- function(df, cloneCall = "gene+nt", chain = "both", scale=FALSE, group = NULL, 
                    exportTable = FALSE) {
    if (length(group) > 1) { stop("Only one item in the group variable can 
                                    be listed.") }
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    if (chain != "both") {
      for(i in seq_along(df)) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
    }
    if (!is.null(group)) {
        x <- group
        labs <- group
        Con.df <- data.frame(matrix(NA, length(df), 4))
        colnames(Con.df) <- c("contigs","values", "total", group)
        for (i in seq_along(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall])
            location <- which(colnames(df[[i]]) == group)
            Con.df[i,4] <- df[[i]][1,location] }
        col <- length(unique(Con.df[,group]))
        if (scale == TRUE) { y <- "scaled"
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"

        } else { y <- "contigs"
            x <- group
            ylab <- "Unique Clonotypes"}
    } else {
        x <- "values"
        labs <- "Samples"
        Con.df <- data.frame(matrix(NA, length(df), 3))
        colnames(Con.df) <- c("contigs","values", "total")
        for (i in seq_along(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall]) }
        col <- length(unique(Con.df$values))
        if (scale == TRUE) { y <- "scaled"
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"
        } else { y <- "contigs"
            ylab <- "Unique Clonotypes" } }
    if (exportTable == TRUE) { return(Con.df) }
    plot <- ggplot(aes(x=Con.df[,x], y=Con.df[,y],
            fill=as.factor(Con.df[,x])), data = Con.df) +
        stat_summary(geom = "errorbar", fun.data = mean_se, 
            position = "dodge", width=.5) + labs(fill = labs) +
        stat_summary(fun=mean, geom="bar", color="black", lwd=0.25)+
        theme_classic() + xlab("Samples") + ylab(ylab) +
        scale_fill_manual(values = colorblind_vector(col))
    return(plot) }


#' Demonstrate the relative abundance of clonotypes by group or sample.
#'
#' This function takes the output of combineTCR(), combineBCR(), or 
#' expression2List() and displays the number of clonotypes at specific 
#' frequencies by sample or group. Visualization can either be a line 
#' graph using calculated numbers or if scale = TRUE, the output will 
#' be a density plot. Multiple sequencing runs can be group together 
#' using the group parameter. If a matrix output for the data is 
#' preferred, set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' abundanceContig(combined, cloneCall = "gene", scale = FALSE)
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List().
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param group The column header for which you would like to analyze the data.
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param exportTable Returns the data frame used for forming the graph
#' to the visualization.
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the total or relative abundance of clonotypes 
#' across quanta
abundanceContig <- function(df, cloneCall = "gene+nt", chain = "both", scale=FALSE, 
                        group = NULL, exportTable = FALSE) {
    Con.df <- NULL
    xlab <- "Abundance"
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    names <- names(df)
    if (!is.null(group)) {
        for (i in seq_along(df)) {
            if (chain != "both") {
              df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
            }
            data1 <- parseContigs(df, i, names, cloneCall)
            label <- df[[i]][1,group]
            data1[,paste(group)] <- label
            Con.df<- rbind.data.frame(Con.df, data1) }
        Con.df <- data.frame(Con.df)
        col <- length(unique(Con.df[,group]))
        fill <- group
        if (scale == TRUE) { ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, fill=Con.df[,group])) +
                geom_density(aes(y=..scaled..), alpha=0.5, 
                    lwd=0.25, color="black", bw=0.5)  +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else { ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                    color = Con.df[,group])) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)}
    } else{
        for (i in seq_along(df)) {
          if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
          }
            data1 <- parseContigs(df, i, names, cloneCall)
            Con.df<- rbind.data.frame(Con.df, data1) }
        col <- length(unique(Con.df$values))
        fill <- "Samples"
        if (scale == TRUE) { ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, 
                    color="black", bw=0.5) +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else { ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group = values, 
                    color = values)) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)} }
    if (exportTable == TRUE) { return(Con.df)}
    plot <- plot + scale_x_log10() + ylab(ylab) + xlab(xlab) +
        theme_classic()
return(plot) }

#' Demonstrate the distribution of lengths filtered contigs.
#'
#' This function takes the output of combineTCR(), combineBCR(), or
#' expression2List() and displays either the nucleotide (nt) or amino 
#' acid (aa) sequence length. The sequence length visualized can be 
#' selected using the chains parameter, either the combined clonotype 
#' (both chains) or across all single chains. Visualization can either 
#' be a histogram or if scale = TRUE, the output will be a density plot. 
#' Multiple sequencing runs can be group together using the 
#' group parameter. If a matrix output for the data is preferred, set 
#' exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' lengthContig(combined, cloneCall="aa", chain = "both")
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List()
#' @param cloneCall How to call the clonotype - CDR3 nucleotide (nt), 
#' CDR3 amino acid (aa).
#' @param group The group header for which you would like to analyze 
#' the data.
#' @param scale Converts the graphs into density plots in order to show 
#' relative distributions.
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param exportTable Returns the data frame used for forming the graph.
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot
#' @export
#' @return ggplot of the discrete or relative length distributions of 
#' clonotype sequences
lengthContig <- function(df, cloneCall = "aa", group = NULL, scale = FALSE, 
                    chain = "both", exportTable = FALSE) {
    if(cloneCall == "nt") { cloneCall <- "CTnt"
        ylab <- "CDR3 (NT)"
    } else if (cloneCall == "aa") { cloneCall <- "CTaa" 
        ylab <- "CDR3 (AA)"
    } else { stop("Please make a selection of the type of
                CDR3 sequence to analyze by using `cloneCall`") }
    cells <- df[[1]][1,"cellType"]
    c1 <- cellT(cells)[[1]]
    c2 <- cellT(cells)[[2]]
    xlab <- "Length"
    Con.df <- NULL
   
    Con.df <- lengthDF(df, cloneCall, chain, group, c1, c2)
    names <- names(df)
    if (!is.null(group)) { 
        fill = group
        col <- length(unique(Con.df[,group]))
        if (scale == TRUE) { yplus <- "Percent of "
            plot <- ggplot(Con.df, aes(fill=Con.df[,group],
                length,(..count..)/sum(..count..)*100)) + 
                geom_density(aes(y=..scaled..),alpha=.5,lwd=.25,color="black")
        } else { yplus <- "Number of "
            plot<-ggplot(Con.df,aes(as.factor(length),fill=Con.df[,group]))+
                geom_bar(position = position_dodge2(preserve = "single"), 
                color="black", lwd=0.25, width=0.9)  +
                scale_x_discrete(breaks = round(seq(min(Con.df$length), 
                max(Con.df$length), by = 5),10)) }
    } else if (is.null(group)){ 
        fill <- "Samples"
        col <- length(unique(Con.df$values))
        if (scale == TRUE) { yplus <- "Percent of "
            plot <- ggplot(Con.df, aes(length, (..count..)/sum(..count..)*100, 
                fill=values)) + geom_density(aes(y=..scaled..), alpha=0.5, 
                lwd=0.25, color="black")
    }  else { yplus <- "Number of "
        plot <- ggplot(Con.df, aes(as.factor(length), fill=values)) +
            geom_bar(position = position_dodge2(preserve = "single"), 
            color="black", lwd=0.25) +
            scale_x_discrete(breaks = round(seq(min(Con.df$length), 
            max(Con.df$length), by = 5),10))} }
    plot <- plot + scale_fill_manual(values = colorblind_vector(col)) +
        labs(fill = fill) + ylab(paste(yplus, ylab, sep="")) +
        xlab(xlab) + theme_classic()
    if (exportTable == TRUE) { return(Con.df) }
    return(plot)}

#' Demonstrate the difference in clonal proportion between clonotypes
#'
#' This function produces an alluvial or area graph of the proportion of 
#' the indicated clonotypes for all or selected samples. Clonotypes can be 
#' selected using the clonotypes parameter with the specific sequence of 
#' interest or using the number parameter with the top n clonotypes by 
#' proportion to be visualized. If multiple clonotypes have the same proportion 
#' and are within the selection by the number parameter, all the clonotypes 
#' will be visualized. In this instance, if less clonotypes are desired, 
#' reduce the number parameter.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' compareClonotypes(combined, numbers = 10, 
#' samples = c("PX_P", "PX_T"), cloneCall="aa")
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List()
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param samples The specific samples to isolate for visualization.
#' @param clonotypes The specific sequences of interest.
#' @param numbers The top number clonotype sequences.
#' @param graph The type of graph produced, either "alluvial" or "area".
#' @param exportTable Returns the data frame used for forming the graph.
#' @import ggplot2
#'
#' @export
#' @return ggplot of the proportion of total sequencing read of 
#' selecting clonotypes
compareClonotypes <- function(df, cloneCall = "gene+nt", chain = "both", samples = NULL, 
                        clonotypes = NULL, numbers = NULL, graph = "alluvial",
                        exportTable = FALSE){
    cloneCall <- theCall(cloneCall)
    df <- checkBlanks(df, cloneCall)
    if (!is.null(numbers) & !is.null(clonotypes)) {
        stop("Make sure your inputs are either numbers or clonotype sequences.")
    }
    Con.df <- NULL
    for (i in seq_along(df)) {
      if (chain != "both") {
        df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
      }
        tbl <- as.data.frame(table(df[[i]][,cloneCall]))
        tbl[,2] <- tbl[,2]/sum(tbl[,2])
        colnames(tbl) <- c("Clonotypes", "Proportion")
        tbl$Sample <- names(df[i])
        Con.df <- rbind.data.frame(Con.df, tbl)
    }
    if (!is.null(samples)) {
        Con.df <- Con.df[Con.df$Sample %in% samples,] }
    if (!is.null(clonotypes)) {
        Con.df <- Con.df[Con.df$Clonotypes %in% clonotypes,] }
    if (!is.null(numbers)) {
        top <- Con.df %>% top_n(n = numbers, wt = Proportion)
        Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,] }
    if (nrow(Con.df) < length(unique(Con.df$Sample))) {
        stop("Reasses the filtering strategies here, there is not 
            enough clonotypes to examine.") }
    if (exportTable == TRUE) { return(Con.df)}
    
    plot <- ggplot(Con.df, aes(x = Sample, fill = Clonotypes, group = Clonotypes,
                    stratum = Clonotypes, alluvium = Clonotypes, 
                    y = Proportion, label = Clonotypes)) +
                theme_classic() +
                theme(axis.title.x = element_blank())
    if (graph == "alluvial") {
        plot = plot +  geom_stratum() + geom_flow(stat = "alluvium")
    } else if (graph == "area") {
        plot = plot +
            geom_area(aes(group = Clonotypes), color = "black") }
    return(plot)
}

#' Scatter plot comparing the expansion of two samples
#'
#' This function produces a scatter plot directly comparing the specific clonotypes
#' between two samples. The clonotypes will be categorized by counts into singlets or multiplets, 
#' either exclusive or shared between the selected samples. Visualization inspired 
#' by the work of \href{https://pubmed.ncbi.nlm.nih.gov/32103181/}{Wu, T, et al 2020}. 
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' scatterClonotype(combined, x.axis = "PY_P", y.axis = "PY_T",
#' graph = "proportion")
#' 
#' @param df The product of combineTCR(), combineBCR(), or expression2List()
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param x.axis name of the list element to appear on the x.axis
#' @param y.axis name of the list element to appear on the y.axis
#' @param dot.size either total or the name of the list element to 
#' use for size of dots
#' @param graph graph either proportion or raw clonotype count
#' @param exportTable Returns the data frame used for forming the graph.
#' 
#' @import ggplot2
#' 
#' @export
#' @return ggplot of the relative clonotype numbers

scatterClonotype <- function(df, cloneCall ="gene+nt", 
                             x.axis = NULL, y.axis = NULL,
                             chain = "both",
                             dot.size = "total", 
                             graph = "proportion", 
                             exportTable = FALSE) {
  `%!in%` = Negate(`%in%`)
  cloneCall <- theCall(cloneCall)
  axes <- which(names(df) %in% c(x.axis, y.axis, dot.size))
  if (chain != "both") {
    for (i in axes) {
      df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
    }
  }
  x.df <- as.data.frame(table(df[[x.axis]][,cloneCall]))
  colnames(x.df)[2] <- x.axis
  y.df <- as.data.frame(table(df[[y.axis]][,cloneCall]))
  colnames(y.df)[2] <- y.axis
  combined.df <- merge(x.df, y.df, by = "Var1", all = TRUE)
  if (dot.size != "total") {
    if (dot.size %!in% colnames(combined.df)) {
      size.df <- as.data.frame(table(df[[dot.size]][,cloneCall]))
      colnames(size.df)[2] <- dot.size
      combined.df <- merge(combined.df, size.df, by = "Var1", all = TRUE) }
    combined.df[is.na(combined.df)] <- 0
    combined.df[,paste0("size", ".fraction")] <- combined.df[,dot.size]/sum(combined.df[,dot.size])
    labeling <- unique(c(x.axis, y.axis, dot.size))
  } else {
    combined.df[is.na(combined.df)] <- 0
    labeling <- unique(c(x.axis, y.axis)) }
  combined.df[,"class"] <- NA
  combined.df[,"sum"] <- rowSums(combined.df[,labeling])
  for (i in seq_along(labeling)) {
    if (length(labeling) > 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & rowSums(combined.df[,labeling[which(labeling != labeling[i])]]) == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
    } else if (length(labeling) == 2) {
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] == 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".singlet"), combined.df[,"class"])
      combined.df[,"class"] <- ifelse(combined.df[,labeling[i]] > 1 & combined.df[,labeling[which(labeling != labeling[i])]] == 0, 
                                      paste0(labeling[i], ".multiplet"), combined.df[,"class"])
  }}
  combined.df[,"class"] <- ifelse(combined.df[,y.axis] >= 1 & combined.df[,x.axis] >= 1, paste0("dual.expanded"), combined.df[,"class"])
  combined.df[,paste0(x.axis, ".fraction")] <- combined.df[,x.axis]/sum(combined.df[,x.axis])
  combined.df[,paste0(y.axis, ".fraction")] <- combined.df[,y.axis]/sum(combined.df[,y.axis])
  if (graph == "proportion") {
    x <- combined.df[,paste0(x.axis, ".fraction")]
    y <- combined.df[,paste0(y.axis, ".fraction")]
  } else if (graph == "count") {
    x <- combined.df[,x.axis]
    y <- combined.df[,y.axis] }
  if (dot.size != "total") {
    size <- combined.df[,dot.size]
  } else { size <- combined.df[,"sum"] }
  if (exportTable == TRUE) { return(combined.df) }
  plot <- ggplot(combined.df, aes(x=x, y = y, color = class)) + 
    theme_classic() + 
    scale_color_manual(values = colorblind_vector(length(unique(combined.df$class)))) + 
    xlab(x.axis) + ylab(y.axis) + labs(size = "Total n")
  if (graph == "proportion") {
    plot <- plot + geom_abline(slope = 1, intercept = 0, alpha = 0.4, lty=2)  + 
      scale_y_sqrt() + scale_x_sqrt() 
  } else if (graph == "count") {
    plot <- plot + ylim(0, max(x,y)) + xlim(0, max(x,y)) 
  }
  plot <- plot + geom_jitter(aes(size = size))
  return(plot)  
}

#' Hierarchical clustering of clonotypes on clonotype size and 
#' Jensen-Shannon divergence
#'
#' This function produces a hierarchical clustering of clonotypes by sample 
#' using the Jensen-Shannon distance and discrete gamma-GPD spliced threshold 
#' model in 
#' \href{https://bioconductor.org/packages/devel/bioc/html/powerTCR.html}{powerTCR R package}.
#' Please read and cite PMID: 30485278 if using the function for analyses. 
#' 
#' 
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonesizeDistribution(combined, cloneCall = "gene+nt", method="ward.D2")
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List().
#' @param cloneCall How to call the clonotype - VDJC gene (gene),
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (gene+nt).
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "both", "TRA", "TRG", "IGH", "IGL"
#' @param method The clustering paramater for the dendrogram.
#' @param exportTable Returns the data frame used for forming the graph.
#' @importFrom  dplyr bind_rows
#' @importFrom ggplot2 ggplot
#' @importFrom powerTCR fdiscgammagpd get_distances
#' @export
#' @return ggplot dendrogram of the clone size distribution

clonesizeDistribution <- function(df,  cloneCall ="gene+nt", chain = "both",
                            method = "ward.D2", exportTable = FALSE) {
        cloneCall <- theCall(cloneCall)
        data <- bind_rows(df)
        unique_df <- unique(data[,cloneCall])
        Con.df <- data.frame(matrix(NA, length(unique_df), length(df)))
        Con.df <- data.frame(unique_df, Con.df, stringsAsFactors = FALSE)
        colnames(Con.df)[1] <- "clonotype"
        for (i in seq_along(df)) {
          if (chain != "both") {
            df[[i]] <- off.the.chain(df[[i]], chain, cloneCall)
          }
            data <- df[[i]]
            data <- data.frame(table(data[,cloneCall]), 
                        stringsAsFactors = FALSE)
            colnames(data) <- c(cloneCall, "Freq")
            for (y in seq_along(unique_df)){
                    clonotype.y <- Con.df$clonotype[y]
                    location.y <- which(clonotype.y == data[,cloneCall])
                    Con.df[y,i+1] <- data[location.y[1],"Freq"] }
        }
        colnames(Con.df)[2:(length(df)+1)] <- names(df)
        Con.df[is.na(Con.df)] <- 0
        list <- list()
        for (i in seq_along(df)) {
            list[[i]] <- Con.df[,i+1]
            list[[i]] <- suppressWarnings(fdiscgammagpd(list[[i]], useq = 1))
            }
        names(list) <- names(df)
        grid <- 0:10000
        distances <- get_distances(list, grid, modelType="Spliced")
        hclust <- hclust(as.dist(distances), method = method)
        hcd <- as.dendrogram(hclust)
        plot <- plot(hcd)
        if (exportTable == TRUE) { return(distances) }
        return(plot)
}

#This is the basic color palette for the package
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", 
                                        "#C6FDEC", "#7AC5FF", "#0348A6"))

#Making lodes to function in alluvial plots
#' @importFrom ggalluvial to_lodes_form
makingLodes <- function(meta2, color, alpha, facet, set.axes) {
    if (!is.null(color) & !is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id="alluvium",
        axes=set.axes,diffuse=c(as.name(color),as.name(alpha),as.name(facet)))
    } else  if (!is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id="alluvium",
        axes = set.axes, diffuse = c(as.name(color), as.name(alpha)))
    } else if (!is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2,key="x",value="stratum",id ="alluvium",
        axes =set.axes, diffuse = c(as.name(color), as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id="alluvium",axes=set.axes,diffuse=c(as.name(alpha),as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(facet)))
    } else if (!is.null(color) & is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(color)))
    } else if (is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
        id = "alluvium", axes = set.axes, diffuse = c(as.name(alpha)))
    } else { lodes <- to_lodes_form(meta2, key = "x", value = "stratum", 
            id = "alluvium", axes = set.axes)}
    return(lodes) }


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
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' vizGenes(combined, gene = "V", chain = "TRB", plot = "bar", scale = TRUE)
#'
#' @param df The product of combineTCR(), combineBCR(), or expression2List().
#' @param gene Which part of the immune receptor to visualize - V, D, J, C
#' @param chain indicate if both or a specific chain should be used - 
#' e.g. "TRA", "TRG", "IGH", "IGL" (no both option here)
#' @param plot The type of plot to return - heatmap or bar 
#' @param separate Categorical variable which to separate by along y-axis
#' @param order Categorixal variable to organize the x-axis, either "gene" or "variance"
#' @param scale Converts the proportion of total genes 
#' @param exportTable Returns the data frame used for forming the graph.
#' @import ggplot2
#' @importFrom stringr str_split
#' @export
#' @return ggplot bar diagram of vgene counts

vizGenes <- function(df, 
                      gene = "V",
                      chain = "TRA", 
                      plot = "heatmap",  
                      separate = "sample", 
                      order = "gene",
                      scale = TRUE, 
                      exportTable = FALSE) {
    df <- bind_rows(df, .id = "element.names")
    "%!in%" <- Negate("%in%")
    if (separate %!in% colnames(df) | is.null(separate)) {
      seperate <- "element.names"
      }
    if (chain %in% c("TRB", "TRG", "IGH")) {
        gene <- unname(c(V = 1, D = 2, J = 3, C = 4)[gene])
      } else {
        gene <- unname(c(V = 1, J = 2, C = 3)[gene])
      }
    if (ncol(str_split(df[,"CTgene"], "_", simplify = TRUE)) == 1) {
      C1 <- str_split(df[,"CTgene"], "_", simplify = TRUE)[,1] 
      C1 <- str_split(C1, "[.]", simplify = TRUE)[,gene] 
      df$C1 <- C1
      x <- "C1"
    } else {
        C1 <- str_split(df[,"CTgene"], "_", simplify = TRUE)[,1] 
        C1 <- str_split(C1, "[.]", simplify = TRUE)[,gene] 
        C2 <- str_split(df[,"CTgene"], "_", simplify = TRUE)[,2] 
        C2 <- str_split(C2, "[.]", simplify = TRUE)[,gene] 
        df$C1 <- C1
        df$C2 <- C2
        if (chain %in% c("TRA", "TRD", "IGH")) {
            x <- "C1"}
        else if (chain %in% c("TRB", "TRG", "IGL")) {
            x <- "C2"}
    }
    df <- subset(df, !is.na(df[,x])) #remove NA values
    df <- subset(df, df[,x] != "NA") #remove values that are character "NA"
    df <- subset(df, df[,x] != "") #remove rows with non genes
    df <- table(df[,x], df[,separate])
    if (scale == TRUE) {
      for (i in seq_len(ncol(df))) {
        df[,i] <- df[,i]/sum(df[,i])
      }
    }
    df <- as.data.frame(df)
    if (order == "variance") {
      df2 <- df %>%
        group_by(Var1) %>%
        summarise(varcount = var(Freq))
      geneOrder <- order(df2$mcount, decreasing = TRUE)
      varOrder <- order(df2$varcount, decreasing = TRUE)
      df$Var1 <- factor(df$Var1, levels = df$Var1[varOrder])
    }
    if (plot == "bar") {
      plot <- ggplot(df, aes(x=Var1, y = Freq)) + 
          geom_bar(stat = "identity") + 
          theme_classic() + 
          theme(axis.title.x = element_blank(), #remove titles
                axis.title.y = element_blank(), 
                  axis.ticks.x = element_blank(), #removes ticks
                  axis.text.x = element_text(angle = 90, 
                                      vjust = 0.5, hjust=1, size=rel(0.5))) + 
        facet_grid(Var2~.)
    } else if (plot == "heatmap") {
     
      plot <- ggplot(df, aes(x=Var1, y = Var2, fill = Freq)) + 
      geom_tile() + 
      theme_classic() + 
      theme(axis.title.x = element_blank(), #remove titles
            axis.ticks.x = element_blank(), #removes ticks
            axis.text.x = element_text(angle = 90, 
                                       vjust = 0.5, hjust=1, size=rel(0.5)), 
            axis.title.y = element_blank()) + 
        scale_fill_gradientn(colors = rev(colorblind_vector(5)))
    }
    if (exportTable == TRUE) { return(df) }
    return(plot)
}
