#This is the basic color palette for the package
#' @import RColorBrewer
#' @import colorRamps
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

#' Quantify the unique clonotypes in the filtered contigs output from 10x Genomics
#'
#' @param df The product of CombineContig()
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param column The column header for which you would like to analyze the data
#' @param scale Converts the graphs into percentage of unique clonotypes
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @import ggplot2
#' @export
quantContig <- function(df,
                        cloneCall = c("gene", "nt", "aa", "gene+nt"),
                        scale=F,
                        column = NULL,
                        exportTable = F) {
    if (length(column) > 1) {
        stop("Only one item in the column variable can be listed.")
    }
    cloneCall <- theCall(cloneCall)

    if (!is.null(column)) {
        Con.df <- data.frame(matrix(NA, length(df), 4))
        colnames(Con.df) <- c("contigs","values", "total", column)
        for (i in 1:length(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall])
            location <- which(colnames(df[[i]]) == column)
            Con.df[i,4] <- df[[i]][1,location]
        }
        col <- length(unique(Con.df[,column]))
        if (scale == T) {
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"
            y <- "scaled"
            x <- column
            labs <- column
        } else {
            y <- "contigs"
            x <- column
            ylab <- "Unique Clonotypes"
            labs <- column

        }

    } else {
        Con.df <- data.frame(matrix(NA, length(df), 3))
        colnames(Con.df) <- c("contigs","values", "total")
        for (i in 1:length(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,cloneCall]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,cloneCall])
        }
        col <- length(unique(Con.df$values))
        if (scale == T) {
            Con.df$scaled <- Con.df$contigs/Con.df$total*100
            ylab <- "Percent of Unique Clonotype"
            y <- "scaled"
            x <- "values"
            labs <- "Samples"


        } else {
            y <- "contigs"
            x <- "values"
            ylab <- "Unique Clonotypes"
            labs <- "Samples"
        }
    }
    if (exportTable == T) {
        quantContig_output <<- Con.df
    }
    plot <- ggplot(aes(x=Con.df[,x], y=Con.df[,y],fill=as.factor(Con.df[,x])), data = Con.df) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.5) +
        stat_summary(fun.y=mean, geom="bar", color="black", lwd=0.25)+
        theme_classic() +
        xlab("Samples") +
        ylab(ylab) +
        labs(fill = labs) +
        scale_fill_manual(values = colorblind_vector(col))
    return(plot)
}


#' Demonstrate the relative abundance of filtered contigs output from 10x Genomics
#'
#' @param df The product of CombineContig().
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param column The column header for which you would like to analyze the data
#' @param scale Converts the graphs into denisty plots in order to show relative distributions.
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @importFrom ggplot2 ggplot
#' @export
abundanceContig <- function(df,
                            cloneCall = c("gene", "nt", "aa", "gene+nt"),
                            scale=F,
                            column = NULL,
                            exportTable = F) {
    Con.df <- NULL
    xlab <- "Abundance"
    cloneCall <- theCall(cloneCall)
    if (!is.null(column)) {
        for (i in seq_along(df)) {
            names <- names(df)
            data <- df[[i]]
            data1 <- data %>%
                group_by(data[,cloneCall]) %>%
                summarise(Abundance=n())
            colnames(data1)[1] <- cloneCall
            data1$values <- names[i]
            label <- df[[i]][1,column]
            data1[,paste(column)] <- label
            Con.df<- rbind.data.frame(Con.df, data1)
        }
        Con.df <- data.frame(Con.df)
        col <- length(unique(Con.df[,column]))
        fill <- column
        if (scale == T) {
            ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, fill=Con.df[,column])) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black", bw=0.5)  +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else {
            ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group = values, color = Con.df[,column])) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)
        }
    } else{
        for (i in seq_along(df)) {
            names <- names(df)
            data <- df[[i]]
            data1 <- data %>%
                group_by(data[,cloneCall]) %>%
                summarise(Abundance=n())
            colnames(data1)[1] <- cloneCall
            data1$values <- names[i]
            Con.df<- rbind.data.frame(Con.df, data1)
        }
        col <- length(unique(Con.df$values))
        fill <- "Samples"
        if (scale == T) {
            ylab <- "Density of Clonotypes"
            plot <- ggplot(Con.df, aes(Abundance, fill=values)) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black", bw=0.5) +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)

        } else {
            ylab <- "Number of Clonotypes"
            plot <- ggplot(Con.df, aes(x=Abundance, group = values, color = values)) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)
        }
    }
    if (exportTable == T) {
        abundanceContig_output <<- Con.df
    }
    plot <- plot +
        scale_x_log10() +
        ylab(ylab) +
        xlab(xlab) +
        theme_classic()
    return(plot)
}

#' Demonstrate the distribution of lengths filtered contigs output from 10x Genomics
#'
#' @description
#' Examine eithe the nucleotide (nt) or amino acid (aa) sequence length across samples or by vairable (column). Sequences can be visualized as combined values (both chains), or as single chains. If more than 2 chains of the same locus are assigned to the individual barcode, this function will take the first chain.
#'
#' @param df The product of CombineContig()
#' @param cloneCall How to call the clonotype - CDR3 nt or CDR3 aa sequence.
#' @param column The column header for which you would like to analyze the data
#' @param scale Converts the graphs into denisty plots in order to show relative distributions.
#' @param chains Whether to keep clonotypes "combined" or visualize by chain
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot
#' @export
lengthContig <- function(df,
                         cloneCall = c("nt", "aa"),
                         column = NULL,
                         scale = F,
                         chains = c("combined", "single"),
                         exportTable = F) {
    if(cloneCall == "nt") {
        cloneCall <- "CTnt"
        ylab <- "CDR3 (NT)"
    } else if (cloneCall == "aa") {
        cloneCall <- "CTaa"
        ylab <- "CDR3 (AA)"
    }
    else {
        stop("Please make a selection of the type of CDR3 sequence to analyze by using `cloneCall`")
    }
    cells <- df[[1]][1,"cellType"]
    if (cells == "T-AB") {
        c1 <- "TCRA"
        c2 <- "TCRB"
    } else if (cells == "T-GD") {
        c1 <- "TCRG"
        c2 <- "TCRD"
    } else if (cells == "B") {
        c1 <- "IGH"
        c2 <- "IGL"
    }
    xlab <- "Length"
    Con.df <- NULL
    if (!is.null(column)) {
        fill = column
        names <- names(df)
        if (chains == "combined") {
            for (i in seq_along(df)) {
                length <- nchar(df[[i]][,cloneCall])
                val <- df[[i]][,cloneCall]
                cols <- df[[i]][,column]
                data <- data.frame(length, val, cols, names[i])
                data <- na.omit(data)
                colnames(data) <- c("length", "CT", column, "values")
                Con.df<- rbind.data.frame(Con.df, data)
            }
        } else if (chains == "single") {
            for (x in seq_along(df)) {
                strings <- df[[x]][,cloneCall]
                strings <- as.data.frame(str_split(strings, "_", simplify = TRUE), stringsAsFactors = F)
                val1 <- strings[,1]
                for (i in 1:length(val1)) {
                    if (grepl(";", val1[i]) == T) {
                        val1[i] <- str_split(val1, ";", simplify = TRUE)[1]
                    }
                    else {
                        next()
                    }
                }
                chain1 <- nchar(val1)
                cols1 <- df[[x]][,column]
                data1 <- data.frame(chain1, val1, names[x], c1, cols1)
                colnames(data1) <- c("length", "CT", "values", "chain", column)
                val2 <- strings[,2]
                for (i in 1:length(val2)) {
                    if (grepl(";", val2[i]) == T) {
                        val2[i] <- str_split(val2, ";", simplify = TRUE)[1]
                    }
                    else {
                        next()
                    }
                }
                chain2 <- nchar(val2)
                cols2 <- df[[x]][,column]
                data2 <- data.frame(chain2, val2, names[x], c2, cols2)
                colnames(data2) <- c("length", "CT", "values", "chain", column)

                data <- rbind(data1, data2)
                data <- na.omit(data)
                data <- subset(data, CT != "NA")
                data <- subset(data, CT != "")
                Con.df<- rbind.data.frame(Con.df, data)

            }
        }
        col <- length(unique(Con.df[,column]))
        if (scale == T) {
            yplus <- "Percent of "
            plot <- ggplot(Con.df, aes(length, (..count..) / sum(..count..) * 100, fill=Con.df[,column])) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black")
        } else {
            yplus <- "Number of "
            plot <- ggplot(Con.df, aes(as.factor(length), fill=Con.df[,column])) +
                geom_bar(position = position_dodge2(preserve = "single"), color="black", lwd=0.25, width=0.9)  +
                scale_x_discrete(breaks = round(seq(min(Con.df$length), max(Con.df$length), by = 5),10))
        }
    } else if (is.null(column)){
        fill <- "Samples"
        names <- names(df)
        if (chains == "combined") {
            for (i in seq_along(df)) {
                length <- nchar(df[[i]][,cloneCall])
                val <- df[[i]][,cloneCall]
                data <- data.frame(length, val, names[i])
                data <- na.omit(data)
                colnames(data) <- c("length", "CT", "values")
                Con.df<- rbind.data.frame(Con.df, data)
            }
        }
        else if(chains == "single") {
            for (x in seq_along(df)) {
                strings <- df[[x]][,cloneCall]
                strings <- as.data.frame(str_split(strings, "_", simplify = TRUE), stringsAsFactors = F)
                val1 <- strings[,1]
                for (i in 1:length(val1)) {
                    if (grepl(";", val1[i]) == T) {
                        val1[i] <- str_split(val1, ";", simplify = TRUE)[1]
                    }
                    else {
                        next()
                    }
                }
                chain1 <- nchar(val1)
                data1 <- data.frame(chain1, val1, names[x], c1)
                colnames(data1) <- c("length", "CT", "values", "chain")
                val2 <- strings[,2]
                for (i in 1:length(val2)) {
                    if (grepl(";", val2[i]) == T) {
                        val2[i] <- str_split(val2, ";", simplify = TRUE)[1]
                    }
                    else {
                        next()
                    }
                }
                chain2 <- nchar(val2)
                data2 <- data.frame(chain2, val2, names[x], c2)
                colnames(data2) <- c("length", "CT", "values", "chain")

                data <- rbind(data1, data2)
                data <- na.omit(data)
                data <- subset(data, CT != "NA")
                data <- subset(data, CT != "")
                Con.df<- rbind.data.frame(Con.df, data)
            }
        }
    col <- length(unique(Con.df$values))
    if (scale == T) {
        yplus <- "Percent of "
        plot <- ggplot(Con.df, aes(length, (..count..) / sum(..count..) * 100, fill=values)) +
            geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black")
    }  else {
        yplus <- "Number of "
        plot <- ggplot(Con.df, aes(as.factor(length), fill=values)) +
            geom_bar(position = position_dodge2(preserve = "single"), color="black", lwd=0.25) +
            scale_x_discrete(breaks = round(seq(min(Con.df$length), max(Con.df$length), by = 5),10))
    }
    }
    if (chains == "single") {
        plot <- plot + facet_grid(chain~.)
    }
    plot <- plot + scale_fill_manual(values = colorblind_vector(col)) +
        labs(fill = fill) +
        ylab(paste(yplus, ylab, sep="")) +
        xlab(xlab) +
        theme_classic()
    if (exportTable == T) {
        lengthContig_output <<- Con.df
    }
    return(plot)

}

#' Demonstrate the difference in clonal proportion between multiple clonotypes
#'
#' @description
#' Allows for the exaimination of the proportions of selected clonotypes between samples. Specific sequences (clonotypes) can be visualize or users can select the visualization of the top n clonotypes (numbers).
#'
#' @param df The product of CombineContig()
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param samples The specific samples to isolate for visualization
#' @param clonotypes The specific sequences of interest
#' @param numbers The top number clonotype sequences
#' @param graph The type of graph produced, either "alluvial" or "area"
#' @import ggplot2
#'
#' @export
compareClonotypes <- function(df,
                              cloneCall = c("gene", "nt", "aa", "gene+nt"),
                              samples = NULL,
                              clonotypes = NULL,
                              numbers = NULL,
                              graph = c("alluvial", "area")){
    cloneCall <- theCall(cloneCall)
    if (!is.null(numbers) & !is.null(clonotypes)) {
        stop("Make sure your inputs are either numbers or clonotype sequences.")
    }

    Con.df <- NULL
    for (i in seq_along(df)) {
        tbl <- as.data.frame(table(df[[i]][,cloneCall]))
        tbl[,2] <- tbl[,2]/sum(tbl[,2])
        colnames(tbl) <- c("Clonotypes", "Proportion")
        tbl$Sample <- names(df[i])
        Con.df <- rbind.data.frame(Con.df, tbl)
    }

    if (!is.null(samples)) {
        Con.df <- Con.df[Con.df$Sample %in% samples,]
    }
    if (!is.null(clonotypes)) {
        Con.df <- Con.df[Con.df$Sample %in% clonotypes,]
    }
    if (!is.null(numbers)) {
        top <- Con.df %>% top_n(n = numbers, wt = Proportion)
        Con.df <- Con.df[Con.df$Clonotypes %in% top$Clonotypes,]
    }
    if (nrow(Con.df) < length(unique(Con.df$Sample))) {
        stop("Reasses the filtering strategies here, there is not enough clonotypes to examine.")
    }
    plot = ggplot(Con.df, aes(x = Sample, fill = Clonotypes, stratum = Clonotypes,
                              alluvium = Clonotypes, y = Proportion, label = Clonotypes)) +
        theme_classic() +
        theme(axis.title.x = element_blank())

    if (graph == "alluvial") {
        plot = plot +
            geom_flow() +
            geom_stratum()
    } else if (graph == "area") {
        plot = plot +
            geom_area(aes(group = Clonotypes), color = "black")
    }
    return(plot)
}

#' Hierarchical clustering based on JS distance
#'
#' @description
#' Allows for the hierarchical clustering based on Jensen-Shannon distance using the discrete gamma-GPD spliced threshold model in the powerTCR R package. If you are planning on using this function for your analysis, read and cite PMID: 30485278.
#'
#' @param df The product of CombineContig()
#' @param cloneCall How to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param method The clustering paramater for the dendrogram
#' @param exportTable Exports a table of the data into the global environment in addition to the visualization
#' @import dplyr
#' @importFrom ggplot2 ggplot
#' @importFrom powerTCR fdiscgammagpd get_distances
#' @importFrom ggdendro ggdendrogram
#' @export

clonesizeDistribution <- function(df,
                                  cloneCall = c("gene", "nt", "aa", "gene+nt"),
                                  method = c("complete", "ward.D", "ward.D2",
                                             "single", "average", "mcquitty",
                                             "median", "centroid"),
                                  exportTable = F) {
                                                 cloneCall <- theCall(cloneCall)
                                                 data <- bind_rows(df)
                                                 unique_df <- unique(data[,cloneCall])
                                                 Con.df <- data.frame(matrix(NA, length(unique_df), length(df)))
                                                 Con.df <- data.frame(unique_df, Con.df, stringsAsFactors = F)
                                                 colnames(Con.df)[1] <- "clonotype"

                                                 for (i in seq_along(df)) {
                                                     data <- df[[i]]
                                                     data <- data.frame(table(data[,cloneCall]), stringsAsFactors = F)
                                                     colnames(data) <- c(cloneCall, "Freq")
                                                     for (y in 1:length(unique_df)){
                                                         clonotype.y <- Con.df$clonotype[y]
                                                         location.y <- which(clonotype.y == data[,cloneCall])
                                                         Con.df[y,i+1] <- data[location.y[1],"Freq"]
                                                     }
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
                                                 plot <- ggdendrogram(hclust)

                                                 if (exportTable == T) {
                                                     clonesizeDistribution_output <<- distances
                                                 }
                                                 return(plot)
}
