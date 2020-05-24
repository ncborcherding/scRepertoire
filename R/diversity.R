#' Examine the clonal diversity of samples
#'
#' This function calculates traditional measures of diversity - Shannon, 
#' inverse Simpson, Chao1 index, and abundance-based coverage estimators 
#' (ACE) by sample or group. The group paramter can be used to condense 
#' the individual samples. If a matrix output for the data is preferred, 
#' set exportTable = TRUE.
#'
#' @examples
#' #Making combined contig data
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' clonalDiversity(combined, cloneCall = "gene")
#'
#' @param df The product of CombineContig() or expression2List().
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), 
#' CDR3 nucleotide (nt) or CDR3 amino acid (aa), or 
#' CDR3 gene+nucleotide (gene+nt).
#' @param group The column header for which you would like to analyze the data.
#' @param exportTable Exports a table of the data into the global environment 
#' in addition to the visualization
#' @importFrom stringr str_sort
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return ggplot of the diversity of clonotype sequences across list
clonalDiversity <- function(df, cloneCall = c("gene", "nt", "aa", "gene+nt"), 
                            group = c("samples"), exportTable = FALSE) {
    cloneCall <- theCall(cloneCall)
    mat <- NULL
    if (group == "samples") {
        for (i in seq_along(df)) {
            data <- as.data.frame(table(df[[i]][,cloneCall]))
            out <- diversityCall(data)
            mat <- rbind.data.frame(mat, out) }
        colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE")
        rownames(mat) <- names(df)
        mat[,group] <- rownames(mat)
        melt <- melt(mat, id.vars = group)
        plot <- ggplot(melt, aes(x = "", y=value)) +
            geom_jitter(shape=21, size=3, width=0.2, aes(fill=melt[,group]))
        } else {
            for (i in seq_along(df)) {
                data <- as.data.frame(table(df[[i]][,cloneCall]))
                    color <- df[[i]][1,group]
                out <- c(diversityCall(data), color)
                mat <- rbind(mat, out) }
            mat <- as.data.frame(mat)
            colnames(mat) <- c("Shannon", "Inv.Simpson", "Chao", "ACE", group)
            rownames(mat) <- names(df)
            melt <- suppressWarnings(melt(mat, id.vars = group))
            values <- str_sort(as.character(unique(melt[,group])), 
                            numeric = TRUE)
            values2 <- quiet(dput(values))
            melt[,group] <- factor(melt[,group], levels = values2)
            plot <- ggplot(melt, aes(x=melt[,group], y=as.numeric(value))) +
                    geom_jitter(shape=21, size=3, width=0.2, 
                    aes(fill=melt[,group])) }
        col <- length(unique(melt[,group]))
        plot <- plot + ylab("Index Score") + scale_fill_manual(name = group, 
                    values = colorblind_vector(col)) +
                facet_wrap(~variable, scales = "free", ncol = 4) +
                theme_classic() +
                theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank())
        if (exportTable == TRUE) { return(mat) }
        return(plot) }
