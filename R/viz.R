#This is the basic color palette for the package
#' @export
colorblind_vector <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

#df indicates the list of data frames of the contigs
#call is the how to call the clonotype - CDR3 gene, CDR3 nt or CDR3 aa sequence.
#column is the column header for which you would like to analyze the data
#scale converts the graphs into percentage of unique clonotypes

#' @export
quantContig <- function(df,
                        call = c("gene", "nt", "aa"),
                        scale=F,
                        column = NULL) {
    if (length(column) > 1) {
        stop("Only one item in the column variable can be listed.")
    }
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
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
    return(plot)
}

#df indicates the list of data frames of the contigs
#call is the how to call the clonotype - CDR3 genes, CDR3 nt, or CDR3 aa sequence.
#column is the column header for which you would like to analyze the data
#scale converts the graphs into denisty plots in order to show relative distributions.

#' @export
abundanceContig <- function(df,
                            call = c("gene", "nt", "aa"),
                            scale=F,
                            column = NULL) {
    Con.df <- NULL
    xlab <- "Abundance"
    if (call == "gene") {
        call <- "CTgene"
    } else if(call == "nt") {
        call <- "CTnt"
    } else {
        call <- "CTaa"
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
    return(plot)
}

#df indicates the list of data frames of the contigs
#call is the how to call the clonotype - CDR3 nt or CDR3 aa sequence. As of now this only calls the combination of the two loci - producing a bimodal curve for which the smaller is the produced from a single sequence of one of the TCR/Ig loci.
#column is the column header for which you would like to analyze the data
#scale converts the graphs into denisty plots in order to show relative distributions.

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

    return(plot)
}
