#This is the basic color palette for the package
#' @export
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

#' Quantify the unique clonotypes in the filtered contigs output from 10x Genomics
#'
#' @param df The product of CombineContig()
#' @param call is the how to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param column is the column header for which you would like to analyze the data
#' @param scale converts the graphs into percentage of unique clonotypes
#' @example quantContig(combined, call= "gene", scale=T)
#' @export
quantContig <- function(df,
                        call = c("gene", "nt", "aa", "gene+nt"),
                        scale=F,
                        column = NULL) {
    if (length(column) > 1) {
        stop("Only one item in the column variable can be listed.")
    }
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else if (call == "aa") {
        call <- "CTaa"
    } else if (call == "gene+nt") {
        call <- "CTstrict"
    } else {
        stop("Are you sure you made the right call? ", .call = F)
    }
    if (!is.null(column)) {
        Con.df <- data.frame(matrix(NA, length(df), 4))
        colnames(Con.df) <- c("contigs","values", "total", column)
        for (i in 1:length(df)) {
            Con.df[i,1] <- length(unique(df[[i]][,call]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,call])
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
            Con.df[i,1] <- length(unique(df[[i]][,call]))
            Con.df[i,2] <- names(df)[i]
            Con.df[i,3] <- length(df[[i]][,call])
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
    plot <- ggplot2::ggplot(aes(x=Con.df[,x], y=Con.df[,y],fill=as.factor(Con.df[,x])), data = Con.df) +
        stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge", width=.5) +
        stat_summary(fun.y=mean, geom="bar", color="black", lwd=0.25)+
        theme_classic() +
        xlab("Samples") +
        ylab(ylab) +
        labs(fill = labs) +
        scale_fill_manual(values = colorblind_vector(col))
    suppressWarnings(print(plot))
}


#' Demonstrate the relative abundance of filtered contigs output from 10x Genomics
#'
#' @param df The product of CombineContig()
#' @param call is the how to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa, or CDR3+nucleotide
#' @param column is the column header for which you would like to analyze the data
#' @param scale converts the graphs into denisty plots in order to show relative distributions.
#' @example abundanceContig(combined, call= "gene", scale=T)
#'
#' @export
abundanceContig <- function(df,
                            call = c("gene", "nt", "aa", "gene+nt"),
                            scale=F,
                            column = NULL) {
    Con.df <- NULL
    xlab <- "Abundance"
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else if (call == "aa") {
        call <- "CTaa"
    } else if (call == "gene+nt") {
        call <- "CTstrict"
    } else {
        stop("Are you sure you made the right call? ", .call = F)
    }
    if (!is.null(column)) {
        for (i in seq_along(df)) {
            names <- names(df)
            data <- df[[i]]
            data1 <- data %>%
                group_by(data[,call]) %>%
                summarise(Abundance=n())
            colnames(data1)[1] <- call
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
            plot <- ggplot2::ggplot(Con.df, aes(x=Abundance, fill=Con.df[,column])) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black", bw=0.5)  +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)
        } else {
            ylab <- "Number of Clonotypes"
            plot <- ggplot2::ggplot(Con.df, aes(x=Abundance, group = values, color = Con.df[,column])) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)
        }
    } else{
        for (i in seq_along(df)) {
            names <- names(df)
            data <- df[[i]]
            data1 <- data %>%
                group_by(data[,call]) %>%
                summarise(Abundance=n())
            colnames(data1)[1] <- call
            data1$values <- names[i]
            Con.df<- rbind.data.frame(Con.df, data1)
        }
        col <- length(unique(Con.df$values))
        fill <- "Samples"
        if (scale == T) {
            ylab <- "Density of Clonotypes"
            plot <- ggplot2::ggplot(Con.df, aes(Abundance, fill=values)) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black", bw=0.5) +
                scale_fill_manual(values = colorblind_vector(col)) +
                labs(fill = fill)

        } else {
            ylab <- "Number of Clonotypes"
            plot <- ggplot2::ggplot(Con.df, aes(x=Abundance, group = values, color = values)) +
                geom_line(stat="count") +
                scale_color_manual(values = colorblind_vector(col)) +
                labs(color = fill)
        }
    }
    plot <- plot +
        scale_x_log10() +
        ylab(ylab) +
        xlab(xlab) +
        theme_classic()
    suppressWarnings(print(plot))
}

#' Demonstrate the distribution of lengths filtered contigs output from 10x Genomics
#'
#' @param df The product of CombineContig()
#' @param call is the how to call the clonotype - CDR3 nt or CDR3 aa sequence.
#' @param column is the column header for which you would like to analyze the data
#' @param scale converts the graphs into denisty plots in order to show relative distributions.
#' @example lengthContig(combined, call= "aa")
#'
#' @export
lengthContig <- function(df,
                         call = c("nt", "aa"),
                         column = NULL,
                         scale = F) {
    if(call == "nt") {
        call <- "CTnt"
        ylab <- "CDR3 (NT)"
    } else {
        call <- "CTaa"
        ylab <- "CDR3 (AA)"
    }
    xlab <- "Length"
    Con.df <- NULL
    if (!is.null(column)) {
        fill = column
        names <- names(df)
        for (i in seq_along(df)) {
            length <- nchar(df[[i]][,call])
            val <- df[[i]][,call]
            cols <- df[[i]][,column]
            data <- data.frame(length, val, cols, names[i])
            data <- na.omit(data)
            colnames(data) <- c("length", "CT", column, "values")
            Con.df<- rbind.data.frame(Con.df, data)
        }
        col <- length(unique(Con.df[,column]))
        if (scale == T) {
            yplus <- "Percent of "
            plot <- ggplot2::ggplot(Con.df, aes(length, (..count..) / sum(..count..) * 100, fill=Con.df[,column])) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black")
        } else {
            yplus <- "Number of "
            plot <- ggplot2::ggplot(Con.df, aes(as.factor(length), fill=Con.df[,column])) +
                geom_bar(position = position_dodge2(preserve = "single"), color="black", lwd=0.25, width=0.9)  +
                scale_x_discrete(breaks = round(seq(min(Con.df$length), max(Con.df$length), by = 5),10))
        }
    } else{
        fill <- "Samples"
        names <- names(df)
        for (i in seq_along(df)) {
            length <- nchar(df[[i]][,call])
            val <- df[[i]][,call]
            data <- data.frame(length, val, names[i])
            data <- na.omit(data)
            colnames(data) <- c("length", "CT", "values")
            Con.df<- rbind.data.frame(Con.df, data)
        }
        col <- length(unique(Con.df$values))
        if (scale == T) {
            yplus <- "Percent of "
            plot <- ggplot2::ggplot(Con.df, aes(length, (..count..) / sum(..count..) * 100, fill=values)) +
                geom_density(aes(y=..scaled..), alpha=0.5, lwd=0.25, color="black")

        } else {
            yplus <- "Number of "
            plot <- ggplot2::ggplot(Con.df, aes(as.factor(length), fill=values)) +
                geom_bar(position = position_dodge2(preserve = "single"), color="black", lwd=0.25) +
                scale_x_discrete(breaks = round(seq(min(Con.df$length), max(Con.df$length), by = 5),10))
        }
    }
    plot <- plot + scale_fill_manual(values = colorblind_vector(col)) +
        labs(fill = fill) +
        ylab(paste(yplus, ylab, sep="")) +
        xlab(xlab) +
        theme_classic()

    suppressWarnings(print(plot))
}

#' Demonstrate the difference in clonal proportion between multiple clonotypes
#' @param df The product of CombineContig()
#' @param call is the how to call the clonotype - CDR3 nt or CDR3 aa sequence.
#' @param clonotypes the specific sequences of interest
#' @param numbers the top n clonotype sequences
#' @param graph either "alluvial" or "area
#' @example compareContig(combined, call= "aa", numbers=10, graph="alluvial")
#'
#' @export
compareClonotypes <- function(df,
                              call = c("gene", "nt", "aa", "gene+nt"),
                              samples = NULL,
                              clonotypes = NULL,
                              numbers = NULL,
                              graph = c("alluvial", "area")){
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else if (call == "aa") {
        call <- "CTaa"
    } else if (call == "gene+nt") {
        call <- "CTstrict"
    } else {
        stop("Are you sure you made the right call? ", .call = F)
    }
    if (!is.null(numbers) & !is.null(clonotypes)) {
        stop("Make sure your inputs are either numbers or clonotype sequences.")
    }

    Con.df <- NULL
    for (i in seq_along(df)) {
        tbl <- as.data.frame(table(df[[i]][,call]))
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
    suppressWarnings(print(plot))
}
