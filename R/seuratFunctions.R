#' Adding clonotype information to a seurat object
#'
#' @description Add the immune receptor information to the suerat object. By defualt this function also calculates the frequencies of the clonotypes by sequencing run (groupBy = "none"). To change how the frequencies are calculated, select a column header for the groupBy variable.
#' Importantly, before using combineSeurat() ensure the barcodes of the seurat object match the barcodes in the output of the combinedContig() call. Check changeNames() to change the prefix of the seurat object.
#' If the dominant clonotypes have a greater frequency than 500, adjust the cloneTypes variable.
#'
#' @param df The product of CombineTCR() or CombineBCR().
#' @param seurat The seurat object to attach
#' @param cloneCall How to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence, or the combination of genes and nucleotide sequence
#' @param groupBy The column label in the combined contig object in which clonotype frequency will be calculated.
#' @param cloneTypes The bins for the grouping based on frequency
#' @param filterNA Method to subset seurat object of barcodes without clonotype information
#' @import Seurat
#' @export
combineSeurat <- function(df,
                          seurat,
                          cloneCall = c("gene", "nt", "aa", "gene+nt"),
                          groupBy = c("none", "sample", "ID"),
                          cloneTypes = c(Single = 1, Small = 5, Medium = 20, Large = 100, Hyperexpanded = 500),
                          filterNA = F) {
    df <- if(is(df)[1] != "list") list(df) else df
    cloneTypes <- c(None = 0, cloneTypes)
    if (length(df) == 0 | length(seurat) == 0) {
        stop("Make sure you are adding the combined contigs and Seurat object to the combineSeurat() function")
    }
    cloneCall <- theCall(cloneCall)
    Con.df <- NULL
    if (groupBy == "none") {
        for (i in seq_along(df)) {
            data <- data.frame(df[[i]], stringsAsFactors = F)
            data2 <- data[,c("barcode", cloneCall)]
            data2 <- unique(data2)
            data2 <- data2 %>%
                group_by(data2[,cloneCall]) %>%
                summarise(Frequency = n())
            colnames(data2)[1] <- cloneCall
            data2 <- na.omit(data2)
            data <- merge(data, data2, by = cloneCall, all = T)
            Con.df <- rbind.data.frame(Con.df, data)
            data <- NULL

        }
    }
    else if (groupBy != "none") {

        data <- data.frame(bind_rows(df), stringsAsFactors = F)
        data2 <- data[,c("barcode", cloneCall, groupBy)]
        data2 <- unique(data2)
        data2 <- data2 %>%
            group_by(data2[,cloneCall], data2[,groupBy]) %>%
            summarise(Frequency = n())
        data2 <- data.frame(data2, stringsAsFactors = F)
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
    seurat <- AddMetaData(seurat, PreMeta)
    if (filterNA == TRUE) {
        meta <- seurat[[]]
        evalNA <- data.frame(meta[,"cloneType"])

        colnames(evalNA) <- "indicator"
        evalNA <- evalNA %>%
            transmute(indicator = ifelse(is.na(indicator), 0, 1))
        rownames(evalNA) <- rownames(meta)
        seurat <- AddMetaData(seurat, evalNA)
        seurat <- subset(seurat, cloneType != 0)
    }

    return(seurat)
}
#' Highlighting specific clonotypes in Seurat
#'
#' @description Use a specific clonotype sequence to highlight on top of the dimensional reduction in seurat object.
#'
#' @param seurat The seurat object to attach
#' @param cloneCall How to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence or the combination of genes and nucleotide sequence
#' @param sequence The specifc sequence or sequence to highlight
#'
#' @export
highlightClonotypes <- function(seurat,
                                cloneCall = c("gene", "nt", "aa", "gene+nt"),
                                sequence = NULL){
    if (inherits(x=df, what ="Seurat")) {
        stop("Object indicated is not of class 'Seurat', make sure you are using the correct data.")
    }
    if (cloneCall == "gene") {
        cloneCall <- "CTgene"
    } else if(cloneCall == "nt") {
        cloneCall <- "CTnt"
    } else if (cloneCall == "aa") {
        cloneCall <- "CTaa"
    } else if (cloneCall == "gene+nt") {
        cloneCall <- "CTstrict"
    } else {
        stop("Are you sure you made the right cloneCall? ", .cloneCall = F)
    }
    meta <- seurat[[]]
    meta$highlight <- NA
    for(i in seq_along(sequence)) {
        meta$highlight <- ifelse(meta[,cloneCall] == sequence[i], paste("Clonotype", i, sep=""),  meta$highlight)
    }
    names <- rownames(meta)
    meta <- data.frame(meta[,c("highlight")])
    rownames(meta) <- names
    colnames(meta)[1] <- "highlight"
    seurat <- AddMetaData(seurat, meta)

}

#' Exploring interaction of clonotypes by seurat dynamics
#'
#' @description
#' View the proportioal contribution of clonotypes by seurat object meta data. The variable "compare" allows for the selection of the groups, while faceting can seperate the graph.
#'
#' @param seurat The seurat object to attach
#' @param cloneCall How to call the clonotype - based on genes, CDR3 nt, or CDR3 aa sequence or the combination of genes and nucleotide sequence
#' @param compare The column label to visualize
#' @param facet The column label to seperate
#'
#' @import ggfittext
#' @import ggalluvial
#'
#' @export
alluvialClonotypes <- function(seurat,
                          cloneCall = c("gene", "nt", "aa", "gene+nt"),
                          compare = "cluster",
                          facet = NULL) {

    if (inherits(x=df, what ="Seurat")){
        stop("Object indicated is not of class 'Seurat', make sure you are using the correct data.")
    }
    cloneCall <- theCall(cloneCall)
    if (length(compare) > 1) {
        stop("Only one comparison can be made for the sake of sanity, if you'd like to seperate by another variable, use the facet call.")
    }



    meta <- data.frame(seurat[[]], Idents(seurat))
    colnames(meta)[ncol(meta)] <- "cluster"

    if (is.na(meta[,compare])) {
        stop("Make sure you are using the right variable name.")
    }

    plot <- ggplot(meta,
                   aes(axis1 = meta[,compare], axis2 = reorder(meta[,cloneCall], Frequency))) +
        geom_alluvium(aes(fill=meta[,compare])) +
        geom_stratum(fill = "grey", color = "black", lwd=0.1) +
        theme_classic() +
        labs(fill = compare) +
        ylab(compare) +
        geom_fit_text(aes(label = meta[,compare]), stat = "stratum", infer.label = TRUE, reverse = T, min.y=min(table(meta[,compare]))) +
        guides(fill=F) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    if (length(facet) == 1 & length(facet) < 2) {
        plot <- plot + facet_wrap(.~meta[,facet], scales="free_y")
    }
    else if (length(facet) == 0) {
        plot <- plot
    }
    else {
        stop("Only can facet on one item at a time, try multiple calls with different facets.")
    }
    return(plot)
}

#' Change names of the cells in the Seurat object in order to use combineSeurat()
#'
#' @description Change the prefix of the seurat object in order to use combineSeurat().
#'
#' @param seurat The Seurat object
#' @param seuratID The strings/prefixes of the current cell names to replace, i.e. c("PY_Tumor", "PX_Peripheral", "PZ_Tumor")
#' @param newID The strings to replace the Seurat cell names, i.e., c("PY_T", "PX_P", "PZ_T")
#' @export
changeNames <- function(seurat,
                        seuratID = NULL,
                        newID = NULL) {
    if (inherits(x=df, what ="Seurat")) {
        stop("This is to change the cell names of the seurat object.")
    }
    if (length(seuratID) != length(newID)) {
        stop("Make sure the length of the old Seurat prefixes match the length of the prefixes you are changing them with.")
    }
    else {
        x <- rownames(seurat[[]])
        for (i in seq_along(seuratID)) {
            x <- gsub(seuratID[i], newID[i], x)
        }
        seurat <- RenameCells(
            object = seurat,
            new.names = x)

    }
    return(seurat)
}

#' Diversity indices for single-cell RNA-seq
#'
#' @description This function utilizes the Startrac R package derived from PMID: 30479382. Required to run the function, the "type" variable needs to include the difference in where the cells were derived. The output of this function will produce 3 indices: expa (clonal expansion), migra (cross-tissue migration), and trans (state transition). In order to understand the underlying analyses of the outputs please read the manuscript.
#'
#' @param seurat The Seurat object
#' @param type The column header in the meta data that gives the where the cells were derived from, not the patient sample IDs
#' @param sample The column header corresponding to individual samples or patients.
#' @param by Method to subset the indices by either overall (across all samples) or by specific group
#' @importFrom Startrac Startrac.run
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
StartracDiversity <- function(seurat,
                              type = "Type",
                              sample = NULL,
                              by = c("overall")) {
    meta <- data.frame(seurat[[]], Idents(seurat), stringsAsFactors = F)
    colnames(meta)[ncol(meta)] <- "majorCluster"
    meta$clone.status <- ifelse(meta$Frequency > 1, "Yes", "No")
    if (is.null(sample)) {
        stop("Must Add the sample information in order to make the StarTrac calculations")
    } else {
        processed <- data.frame(rownames(meta), meta$CTstrict, meta$clone.status, meta[,sample], meta[,"majorCluster"], meta[,type], stringsAsFactors = F)
        colnames(processed) <- c("Cell_Name", "clone.id", "clone.status", "patient", "majorCluster", "loc")
    }
    processed <- na.omit(processed)
    indices <- suppressWarnings(Startrac.run(processed, proj = "total", verbose = F))
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
            guides(fill=F) +
            theme(axis.title.x = element_blank())

    } else {
        indices <- subset(indices, aid == by)
        plot <- ggplot(melted, aes(x=majorCluster, y=value)) +
            geom_point(aes(fill = by)) +
            facet_grid(variable ~.) +
            theme_classic() +
            ylab("Index Score") +
            guides(fill=F) +
            theme(axis.title.x = element_blank())
    }

    return(plot)

}
