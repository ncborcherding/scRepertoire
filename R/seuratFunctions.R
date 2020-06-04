#' Adding clonotype information to a seurat or SCE object
#'
#' This function adds the immune receptor information to the seurat or 
#' SCE object to the meta data. By defualt this function also calculates 
#' the frequencies of the clonotypes by sequencing run (groupBy = "none"). 
#' To change how the frequencies are calculated, select a column header for 
#' the groupBy variable. Importantly, before using combineExpression() 
#' ensure the barcodes of the seurat or SCE object match the barcodes in the 
#' output of the combinedContig() call. Check changeNames() to change the 
#' prefix of the seurat object. If the dominant clonotypes have a greater 
#' frequency than 500, adjust the cloneTypes variable.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' seurat_example <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/seurat_example.rds"))
#' 
#' #Using combineExpresion()
#' seurat_example <- combineExpression(combined, seurat_example)
#' 
#' @param df The product of CombineTCR() or CombineBCR().
#' @param sc The seurat or SingleCellExperiment (SCE) object to attach
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), 
#' CDR3 nucleotide (nt) CDR3 amino acid (aa), or 
#' CDR3 gene+nucleotide (gene+nt).
#' @param groupBy The column label in the combined contig object in which 
#' clonotype frequency will be calculated.
#' @param cloneTypes The bins for the grouping based on frequency
#' @param filterNA Method to subset seurat object of barcodes without 
#' clonotype information
#' @import Seurat
#' @export
#' @return seurat or SingleCellExperiment object with attached clonotype 
#' information
#' 

combineExpression <- function(df, sc, cloneCall="gene+nt", groupBy="none", 
                        cloneTypes=c(None=0, Single=1, Small=5, Medium=20, 
                        Large=100, Hyperexpanded=500), filterNA = FALSE) {
    df <- scRepertoire:::checkList(df)
    cloneCall <- scRepertoire:::theCall(cloneCall)
    Con.df <- NULL
    meta <- grabMeta(sc)
    cell.names <- rownames(meta)
    if (groupBy == "none") {
        for (i in seq_along(df)) {
            data <- data.frame(df[[i]], stringsAsFactors = FALSE)
            data2 <- unique(data[,c("barcode", cloneCall)])
            data2 <- na.omit(data2[data2[,"barcode"] %in% cell.names,])
            data2 <- data2 %>% group_by(data2[,cloneCall]) %>%
                summarise(Frequency = n())
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
            Con.df <- rbind.data.frame(Con.df, data) }
    } else if (groupBy != "none") {
        data <- data.frame(bind_rows(df), stringsAsFactors = FALSE)
        data2 <- na.omit(unique(data[,c("barcode", cloneCall, groupBy)]))
        data2 <- data2[data2[,"barcode"] %in% cell.names, ]
        data2 <- as.data.frame(data2 %>% group_by(data2[,cloneCall], 
                    data2[,groupBy]) %>% summarise(Frequency = n()))
        colnames(data2)[c(1,2)] <- c(cloneCall, groupBy)
        x <- unique(data[,groupBy])
        for (i in seq_along(x)) {
            sub1 <- subset(data, data[,groupBy] == x[i])
            sub2 <- subset(data2, data2[,groupBy] == x[i])
            merge <- merge(sub1, sub2, by=cloneCall)
            Con.df <- rbind.data.frame(Con.df, merge) } }
    Con.df$cloneType <- NA
    for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <- 
        paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], 
        ' < X <= ', cloneTypes[x], ')') }
    for (i in 2:length(cloneTypes)) { Con.df$cloneType <- 
        ifelse(Con.df$Frequency > cloneTypes[i-1] & Con.df$Frequency 
        <= cloneTypes[i], names(cloneTypes[i]), Con.df$cloneType) }
    PreMeta <- unique(Con.df[,c("barcode", "CTgene", "CTnt", 
                "CTaa", "CTstrict", "Frequency", "cloneType")])
    rownames(PreMeta) <- PreMeta$barcode
    if (inherits(x=sc, what ="Seurat")) { sc <- AddMetaData(sc, PreMeta) 
    } else if (inherits(x=sc, what ="cell_data_set")){
      rownames <- rownames(colData(sc))
      colData(sc) <- 
        merge(colData(sc), PreMeta)[, union(names(colData(sc)), 
                                                 names(PreMeta))]
      rownames(colData(sc)) <- rownames 
    }else{
      rownames <- rownames(sc@metadata[[1]])
      sc@metadata[[1]] <- 
        merge(sc@metadata[[1]], PreMeta)[, union(names(sc@metadata[[1]]), 
                                                 names(PreMeta))]
      rownames(sc@metadata[[1]]) <- rownames 
        }
    if (filterNA == TRUE) { sc <- scRepertoire:::filteringNA(sc) }
    return(sc) }

#' Highlighting specific clonotypes in Seurat
#'
#' Use a specific clonotype sequence to highlight on top of the dimensional 
#' reduction in seurat object.
#'
#' @examples
#' #' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' seurat_example <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/seurat_example.rds"))
#' 
#' #Using combineExpresion()
#' seurat_example <- combineExpression(combined, seurat_example)
#' 
#' #Using highlightClonotype()
#' seurat_example <- highlightClonotypes(seurat_example, cloneCall= "aa", 
#' sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF"))
#' 
#' @param sc The seurat object to attach
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' CDR3 gene+nucleotide (gene+nt).
#' @param sequence The specifc sequence or sequence to highlight
#'
#' @export
#' @return DimPlot with highlighted clonotypes
highlightClonotypes <- function(sc, 
                            cloneCall = c("gene", "nt", "aa", "gene+nt"), 
                            sequence = NULL){
    if (!inherits(x=sc, what ="Seurat")) {
        stop("Object indicated is not of class 'Seurat', make sure you 
            are using the correct data.") }
    cloneCall <- theCall(cloneCall)
    meta <- sc[[]]
    meta$highlight <- NA
    for(i in seq_along(sequence)) {
        meta$highlight <- ifelse(meta[,cloneCall] == sequence[i], 
                            paste("Clonotype", i, sep=""),  meta$highlight) }
    names <- rownames(meta)
    meta <- data.frame(meta[,c("highlight")])
    rownames(meta) <- names
    colnames(meta)[1] <- "highlight"
    sc <- AddMetaData(sc, meta)
}

#' Exploring interaction of clonotypes by seurat or SCE dynamics
#'
#' View the proportional contribution of clonotypes by seurat or SCE object 
#' meta data after combineExpression(). The visualization is based on the 
#' ggalluvial package, which requires the aesthetics to be part of the axes 
#' that are visualized. Therefore, alpha, facet, and color should be part of 
#' the the axes you wish to view or will add an additional stratum/column to 
#' the end of the graph.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' #Getting a sample of a Seurat object
#' seurat_example <- readRDS(url(
#' "https://ncborcherding.github.io/vignettes/seurat_example.rds"))
#' 
#' #Using combineExpresion()
#' seurat_example <- combineExpression(combined, seurat_example)
#' 
#' #Using alluvialClonotypes()
#' alluvialClonotypes(seurat_example, cloneCall = "gene", 
#' y.axes = c("Patient", "cluster"), color = "cluster")
#' 
#' @param sc The seurat or SCE object to visualize after combineExpression(). 
#' For SCE objects, the cluster variable must be in the meta data under 
#' "cluster".
#' @param cloneCall How to call the clonotype - CDR3 gene (gene), 
#' CDR3 nucleotide (nt) or CDR3 amino acid (aa), or 
#' CDR3 gene+nucleotide (gene+nt).
#' @param y.axes The columns that will seperate the proportional 
#' visualizations.
#' @param color The column header or clonotype(s) to be highlighted.
#' @param facet The column label to seperate.
#' @param alpha The column header to have gradieted opacity.
#'
#' @import ggfittext
#' @import ggalluvial
#' @import dplyr
#' @export
#' @return Alluvial ggplot comparing clonotype distribution across 
#' selected parameters.
alluvialClonotypes <- function(sc, 
                                cloneCall = c("gene", "nt", "aa", "gene+nt"), 
                                y.axes = NULL, color = NULL, alpha = NULL, 
                                facet = NULL) {
    checkSingleObject(sc)
    cloneCall <- theCall(cloneCall)
    if (length(y.axes) == 0) {
        stop("Make sure you have selected the variable(s) to visualize") }
    meta <- grabMeta(sc)
    meta$barcodes <- rownames(meta)
    check <- colnames(meta) == color
    if (length(unique(check)) == 1 & unique(check)[1] == FALSE & !is.null(color)) {
        meta <- meta %>% mutate(H.clonotypes = ifelse(meta[,cloneCall] %in% 
            color, "Selected", "Other"))
        color <- "H.clonotypes" }
    y.axes <- unique(c(y.axes, color, alpha, facet))
    set.axes <- seq_along(y.axes)
    meta2 <- meta[,c(y.axes, color, alpha, facet, cloneCall, "barcodes")]
    meta2 <- unique(na.omit(meta2[!duplicated(as.list(meta2))]))
    lodes <- makingLodes(meta2, color, alpha, facet, set.axes) 
    plot <- ggplot(data = lodes, aes(x = x, stratum = stratum, 
                alluvium = alluvium, label = stratum)) +
        geom_stratum() + theme_classic() +
        geom_fit_text(stat = "stratum", infer.label = FALSE, reverse = TRUE) +
        theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    if (is.null(color) & is.null(alpha)) {
        plot <- plot + geom_alluvium()
    } else if (!is.null(color) & is.null(alpha)) {
        plot <- plot+geom_flow(aes(fill = lodes[,color]), 
            stat = "alluvium", lode.guidance = "forward") + labs(fill = color)
    } else if (is.null(color) & !is.null(alpha)) {
        plot <- plot + geom_flow(aes(alpha = lodes[,alpha]), stat = "alluvium",
                lode.guidance = "forward") + labs(alpha = alpha)
    }else {
        plot <- plot+geom_flow(aes(alpha=lodes[,alpha], fill=lodes[,color]),
                        stat = "alluvium", lode.guidance = "forward") + 
                        labs(fill = color, alpha = alpha) }
    if (length(facet) == 1 & length(facet) < 2) {
        plot <- plot + facet_wrap(.~lodes[,facet], scales="free_y")
    } else if (length(facet) == 0) { plot <- plot }
    return(plot)}

