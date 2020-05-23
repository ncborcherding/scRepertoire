#' Adding clonotype information to a seurat or SCE object
#'
#' This function adds the immune receptor information to the seurat or SCE object to the meta data.
#' By defualt this function also calculates the frequencies of the clonotypes by sequencing run (groupBy = "none").
#' To change how the frequencies are calculated, select a column header for the groupBy variable.
#' Importantly, before using combineExpression() ensure the barcodes of the seurat or SCE object match
#' the barcodes in the output of the combinedContig() call. Check changeNames() to change the prefix of the seurat object.
#' If the dominant clonotypes have a greater frequency than 500, adjust the cloneTypes variable.
#'
#' @examples
#' \donttest{
#' combineExpression(seurat, combined, cloneCall = "gene")
#' }
#' @param df The product of CombineTCR() or CombineBCR().
#' @param sc The seurat or SingleCellExperiment (SCE) object to attach
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) or CDR3 amino acid (aa), or
#' CDR3 gene+nucleotide (gene+nt).
#' @param groupBy The column label in the combined contig object in which clonotype frequency will be calculated.
#' @param cloneTypes The bins for the grouping based on frequency
#' @param filterNA Method to subset seurat object of barcodes without clonotype information
#' @import Seurat
#' @export
#' @return seurat or SingleCellExperiment object with attached clonotype information
combineExpression <- function(df,
                          sc,
                          cloneCall = c("gene", "nt", "aa", "gene+nt"),
                          groupBy = c("none", "sample", "ID"),
                          cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500),
                          filterNA = FALSE) {
    df <- if(is(df)[1] != "list") list(df) else df
    cloneTypes <- c(None = 0, cloneTypes)
    if (length(df) == 0 | length(sc) == 0) {
        stop("Make sure you are adding the combined contigs and Seurat/SCE object to the combineExpression() function")
    }
    cloneCall <- theCall(cloneCall)
    Con.df <- NULL
    if (inherits(x=sc, what ="Seurat")) {
        cell.names <- rownames(sc[[]])
    } else if (inherits(x=sc, what ="SummarizedExperiment")){
        cell.names <- rownames(sc@metadata[[1]])
    }
    if (groupBy == "none") {
        for (i in seq_along(df)) {
            data <- data.frame(df[[i]], stringsAsFactors = FALSE)
            data2 <- data[,c("barcode", cloneCall)]
            data2 <- unique(data2)
            data2 <- data2[data2[,"barcode"] %in% cell.names,]
            data2 <- data2 %>%
                group_by(data2[,cloneCall]) %>%
                summarise(Frequency = n())
            colnames(data2)[1] <- cloneCall
            data2 <- na.omit(data2)
            data <- merge(data, data2, by = cloneCall, all = TRUE)
            Con.df <- rbind.data.frame(Con.df, data)
            data <- NULL

        }
    }
    else if (groupBy != "none") {

        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- data[,c("barcode", cloneCall, groupBy)]
        data2 <- unique(data2)
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- data2 %>%
            group_by(data2[,cloneCall], data2[,groupBy]) %>%
            summarise(Frequency = n())
        data2 <- data.frame(data2, stringsAsFactors = FALSE)
        colnames(data2)[1] <- cloneCall
        colnames(data2)[2] <- groupBy
        data2 <- na.omit(data2)
        x <- unique(data[,groupBy])
#        data2<- data2[ , -which(names(data2) %in% c(groupBy))]
        Con.df <- NULL
        for (i in seq_along(x)) {
            sub1 <- subset(data, data[,groupBy] == x[i])
            sub2 <- subset(data2, data2[,groupBy] == x[i])
            merge <- merge(sub1, sub2, by=cloneCall)
            Con.df <- rbind.data.frame(Con.df, merge)
        }
    }
    Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) {
        names(cloneTypes)[x] <- paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], ' < X <= ', cloneTypes[x], ')')
    }
    for (i in 2:length(cloneTypes)) {
        Con.df$cloneType <- ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType)
    }
    PreMeta <- Con.df[,c("barcode", "CTgene", "CTnt", "CTaa", "CTstrict", "Frequency", "cloneType")]
    PreMeta <- unique(PreMeta)
    rownames(PreMeta) <- PreMeta$barcode

    if (inherits(x=sc, what ="Seurat")) {
        sc <- AddMetaData(sc, PreMeta)
    } else if (inherits(x=sc, what ="SummarizedExperiment")){
        rownames <- rownames(sc@metadata[[1]])
        sc@metadata[[1]] <- merge(sc@metadata[[1]], PreMeta)[, union(names(sc@metadata[[1]]), names(PreMeta))]
        rownames(sc@metadata[[1]]) <- rownames
    }
    if (filterNA == TRUE) {
        if (inherits(x=sc, what ="Seurat")) {
            meta <- sc[[]]
        } else if (inherits(x=sc, what ="SummarizedExperiment")){
            meta <- sc@metadata[[1]]
        }

        evalNA <- data.frame(meta[,"cloneType"])

        colnames(evalNA) <- "indicator"
        evalNA <- evalNA %>%
            transmute(indicator = ifelse(is.na(indicator), 0, 1))
        rownames(evalNA) <- rownames(meta)
        sc <- AddMetaData(sc, evalNA)
        sc <- subset(sc, cloneType != 0)
    }

    return(sc)
}
#' Highlighting specific clonotypes in Seurat
#'
#' Use a specific clonotype sequence to highlight on top of the dimensional reduction in seurat object.
#'
#' @examples
#' \donttest{
#' highlightClonotype(seurat, cloneCall = "gene", sequence = "CAVNGGSQGNLIF_CSAEREDTDTQYF")
#' }
#' @param sc The seurat object to attach
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) or CDR3 amino acid (aa), or
#' CDR3 gene+nucleotide (gene+nt).
#' @param sequence The specifc sequence or sequence to highlight
#'
#' @export
#' @return DimPlot with highlighted clonotypes
highlightClonotypes <- function(sc,
                                cloneCall = c("gene", "nt", "aa", "gene+nt"),
                                sequence = NULL){
    if (!inherits(x=sc, what ="Seurat")) {
        stop("Object indicated is not of class 'Seurat', make sure you are using the correct data.")
    }
    cloneCall <- theCall(cloneCall)

    meta <- sc[[]]
    meta$highlight <- NA
    for(i in seq_along(sequence)) {
        meta$highlight <- ifelse(meta[,cloneCall] == sequence[i], paste("Clonotype", i, sep=""),  meta$highlight)
    }
    names <- rownames(meta)
    meta <- data.frame(meta[,c("highlight")])
    rownames(meta) <- names
    colnames(meta)[1] <- "highlight"
    sc <- AddMetaData(sc, meta)

}

#' Exploring interaction of clonotypes by seurat or SCE dynamics
#'
#' View the proportional contribution of clonotypes by seurat or SCE object meta data after combineExpression().
#' The visualization is based on the ggalluvial package, which requires the aesthetics to be part of the axes that are
#' visualized. Therefore, alpha, facet, and color should be part of the the axes you wish to view or will add an additional
#' stratum/column to the end of the graph.
#'
#' @examples
#' \donttest{
#' alluvialClonotypes(seurat, cloneCall = "gene", y.axes = c("Patient", "cluster"), color = "cluster")
#' }
#' @param sc The seurat or SCE object to visualize after combineExpression(). For SCE objects, the cluster variable
#' must be in the meta data under "cluster".
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), CDR3 nucleotide (nt) or CDR3 amino acid (aa), or
#' CDR3 gene+nucleotide (gene+nt).
#' @param y.axes The columns that will seperate the proportional visualizations.
#' @param color The column header or clonotype(s) to be highlighted.
#' @param facet The column label to seperate.
#' @param alpha The column header to have gradieted opacity.
#'
#' @import ggfittext
#' @import ggalluvial
#' @import dplyr
#'
#' @export
#' @return Alluvial ggplot comparing clonotype distribution across selected parameters.
alluvialClonotypes <- function(sc,
                               cloneCall = c("gene", "nt", "aa", "gene+nt"),
                               y.axes = NULL,
                               color = NULL,
                               alpha = NULL,
                               facet = NULL) {

    if (!inherits(x=sc, what ="Seurat") | inherits(x=sc, what ="SummarizedExperiment")){
        stop("Object indicated is not of class 'Seurat' or 'SummarizedExperiment', make sure you are using the correct data.")
    }
    cloneCall <- theCall(cloneCall)
    if (length(y.axes) == 0) {
        stop("Make sure you have selected the variable(s) to visualize")
    }

    if (inherits(x=sc, what ="Seurat")) {
        meta <- data.frame(sc[[]], Idents(sc))
        colnames(meta)[ncol(meta)] <- "cluster"
    } else {
        meta <- sc@metadata[[1]]
    }
    meta$barcodes <- rownames(meta)


    check <- colnames(meta) == color
    if (length(unique(check)) == 1 & unique(check) == FALSE & !is.null(color)) {
        meta <- meta %>%
            mutate(H.clonotypes = ifelse(meta[,cloneCall] %in% color, "Selected", "Other"))
        color <- "H.clonotypes"
    }
    y.axes <- unique(c(y.axes, color, alpha, facet))

    meta2 <- meta[,c(y.axes, color, alpha, facet, cloneCall, "barcodes")]
    meta2 <- meta2[!duplicated(as.list(meta2))]
    meta2 <- na.omit(meta2)
    meta2 <- unique(meta2)


    if (!is.null(color) & !is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(color), as.name(alpha), as.name(facet)))
    } else  if (!is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(color), as.name(alpha)))
    } else if (!is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(color), as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(alpha), as.name(facet)))
    } else if (is.null(color) & is.null(alpha) & !is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(facet)))
    } else if (!is.null(color) & is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(color)))
    } else if (is.null(color) & !is.null(alpha) & is.null(facet)) {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes), diffuse = c(as.name(alpha)))
    } else {
        lodes <- to_lodes_form(meta2, key = "x", value = "stratum", id = "alluvium",
                               axes = seq_along(y.axes))
    }

    plot <- ggplot(data = lodes, aes(x = x, stratum = stratum, alluvium = alluvium, label = stratum)) +
        geom_stratum() +
        geom_fit_text(stat = "stratum",  infer.label = FALSE, reverse = TRUE) +
        theme_classic() +
        theme(axis.title.x = element_blank(),
              axis.ticks.x = element_blank())
    if (is.null(color) & is.null(alpha)) {
        plot <- plot + geom_alluvium()
    } else if (!is.null(color) & is.null(alpha)) {
        plot <- plot +
            geom_flow(aes(fill = lodes[,color]), stat = "alluvium",
                      lode.guidance = "forward") + labs(fill = color)
    } else if (is.null(color) & !is.null(alpha)) {
        plot <- plot +
            geom_flow(aes(alpha = lodes[,alpha]), stat = "alluvium",
                      lode.guidance = "forward") + labs(alpha = alpha)
    }else {
        plot <- plot +
            geom_flow(aes(alpha = lodes[,alpha], fill = lodes[,color]), stat = "alluvium",
                      lode.guidance = "forward") + labs(fill = color, alpha = alpha)
    }

    if (length(facet) == 1 & length(facet) < 2) {
        plot <- plot + facet_wrap(.~lodes[,facet], scales="free_y")
    } else if (length(facet) == 0) {
        plot <- plot
    }
    return(plot)
}

#' Diversity indices for single-cell RNA-seq
#'
#' This function utilizes the Startrac R package derived from PMID: 30479382. Required to run the function,
#' the "type" variable needs to include the difference in where the cells were derived. The output of this function will
#' produce 3 indices: expa (clonal expansion), migra (cross-tissue migration), and trans (state transition). In order to
#' understand the underlying analyses of the outputs please read the manuscript.
#' 
#' @examples
#' \donttest{
#' StartracDiversity(seurat, type = "Type", sample = "Patient", by = "overall")
#' }
#' @param sc The Seurat or SCE object. For SCE objects, the cluster variable must be in the meta data under "cluster".
#' @param type The column header in the meta data that gives the where the cells were derived from, not the patient
#' sample IDs.
#' @param sample The column header corresponding to individual samples or patients.
#' @param by Method to subset the indices by either overall (across all samples) or by specific group.
#' @importFrom Startrac Startrac.run
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
#' @return ggplot of clonal metrics calculatrd by the Startrac R package.
StartracDiversity <- function(sc,
                              type = "Type",
                              sample = NULL,
                              by = c("overall")) {
    if (inherits(x=sc, what ="Seurat")) {
        meta <- data.frame(sc[[]], Idents(sc), stringsAsFactors = FALSE)
        colnames(meta)[ncol(meta)] <- "majorCluster"
    }
    else if (inherits(x=sc, what ="SummarizedExperiment")) {
        meta <- sc["meta.data"][[1]]
        colnames(meta)["cluster"] <- "majorCluster"
    }
    meta$clone.status <- ifelse(meta$Frequency > 1, "Yes", "No")
    if (is.null(sample)) {
        stop("Must Add the sample information in order to make the StarTrac calculations")
    } else {
        processed <- data.frame(rownames(meta), meta$CTstrict, meta$clone.status, meta[,sample], meta[,"majorCluster"], meta[,type], stringsAsFactors = FALSE)
        colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", "majorCluster", "loc")
    }
    processed <- na.omit(processed)
    indices <- suppressWarnings(Startrac.run(processed, proj = "total", verbose = FALSE))
    indices <- data.frame(indices@cluster.data)
    if (by == "overall") {
        indices <- subset(indices, aid != "total")

        melted <- melt(indices)
        values <- as.character(unique(melted$majorCluster))
        values2 <- quiet(dput(values))
        melted$majorCluster <- factor(melted$majorCluster, levels = values2)

        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_boxplot(aes(fill = majorCluster), outlier.alpha = 0) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill=FALSE) +
            theme(axis.title.x = element_blank())

    } else {
        indices <- subset(indices, aid == by)
        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_point(aes(fill = by)) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill=FALSE) +
            theme(axis.title.x = element_blank())
    }

    return(plot)

}
