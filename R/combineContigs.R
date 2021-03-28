# Adding Global Variables
# data('v_gene','j_gene', 'c_gene', 'd_gene')
utils::globalVariables(c("v_gene", "j_gene", "c_gene", "d_gene", "chain"))

heavy_lines <- c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")
light_lines <- c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")
l_lines <- c("IGLct", "cdr3", "cdr3_nt", "v_gene")
k_lines <- c("IGKct", "cdr3", "cdr3_nt", "v_gene")
h_lines <- c("IGHct", "cdr3", "cdr3_nt", "v_gene")
tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
data1_lines <- c("TCR1", "cdr3", "cdr3_nt")
data2_lines <- c("TCR2", "cdr3", "cdr3_nt")
CT_lines <- c("CTgene", "CTnt", "CTaa", "CTstrict", "cellType")

utils::globalVariables(c("heavy_lines", "light_lines", "l_lines", "k_lines", 
            "h_lines", "tcr1_lines", "tcr2_lines", "data1_lines", 
            "data2_lines", "CT_lines"))

#' Combining the list of T Cell Receptor contigs
#'
#' This function consolidates a list of TCR sequencing results to the level of 
#' the individual cell barcodes. Using the samples and ID parameters, the 
#' function will add the strings as prefixes to prevent issues with repeated 
#' barcodes. The resulting new barcodes will need to match the seurat or SCE 
#' object in order to use, \code{\link{combineExpression}}. Several 
#' levels of filtering exist - remove or filterMulti are parameters that 
#' control how  the function  deals with barcodes with multiple chains 
#' recovered.
#' @examples
#' combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' @param df List of filtered contig annotations from 10x Genomics.
#' @param samples The labels of samples.
#' @param ID The additional sample labeling option.
#' @param cells The type of T cell - T cell-AB or T cell-GD
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 2 
#' corresponding chains with the highest expression for a single barcode.
#' @import dplyr
#' @export
#' @return List of clonotypes for individual cell barcodes
combineTCR <- function(df, samples = NULL, ID = NULL, 
                cells = c("T-AB", "T-GD"), removeNA = FALSE, 
                removeMulti = FALSE, filterMulti = FALSE) {
    df <- checkList(df)
    out <- NULL
    final <- NULL
    checkContigBarcodes(df, samples, ID)
    chain1 <- cellT(cells)[[1]]
    chain2 <- cellT(cells)[[2]]
    cellType <- cellT(cells)[[3]]
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain != "Multi")
        df[[i]] <- subset(df[[i]], chain %in% c(chain1, chain2))
        df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        df[[i]]$sample <- samples[i]
        df[[i]]$ID <- ID[i]
        if (filterMulti == TRUE) { df[[i]] <- filteringMulti(df[[i]]) }
        if (nrow(df[[i]]) == 0) { stop("There are 0 contigs 
                after filtering for celltype.", call. = FALSE) }}
    out <- modifyBarcodes(df, samples, ID)
    for (i in seq_along(out)) { 
        data2 <- out[[i]]
        data2 <- makeGenes(cellType, data2, chain1, chain2)
        unique_df <- unique(data2$barcode)
        Con.df <- data.frame(matrix(NA, length(unique_df), 7))
        colnames(Con.df) <- c("barcode",tcr1_lines, tcr2_lines)
        Con.df$barcode <- unique_df
        Con.df <- parseTCR(Con.df, unique_df, data2)
        Con.df <- assignCT(cellType, Con.df)
        Con.df$cellType <- cells 
        Con.df[Con.df == "NA_NA" | Con.df == "NA_NA_NA_NA"] <- NA 
        data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))], 
            Con.df, by = "barcode")
        data3<-data3[,c("barcode","sample","ID",tcr1_lines,tcr2_lines,
            CT_lines)]
        final[[i]] <- data3 }
    names <- NULL
    for (i in seq_along(samples)) { c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names
    for (i in seq_along(final)){
        final[[i]]<-final[[i]][!duplicated(final[[i]]$barcode),]}
    if (removeNA == TRUE) { final <- removingNA(final)}
    if (removeMulti == TRUE) { final <- removingMulti(final) }
    return(final) }

#' Combining the list of B Cell Receptor contigs
#'
#' This function consolidates a list of BCR sequencing results to the level 
#' of the individual cell barcodes. Using the samples and ID parameters, 
#' the function will add the strings as prefixes to prevent issues with 
#' repeated barcodes. The resulting new barcodes will need to match the 
#' seurat or SCE object in order to use, \code{\link{combineExpression}}. 
#' Unlike combineTCR(), combineBCR() produces a column CTstrict of an 
#' index of nucleotide sequence and the corresponding v-gene. This 
#' index automatically calculates the Hammings distance between 
#' sequences of the same length and will index sequences with > 0.85 
#' normalized Hammings distance with the same ID. If nucleotide 
#' sequences meet the threshold, ":HD" will be added to the CTstrict 
#' column string.
#'
#' @examples
#' #Data derived from the 10x Genomics intratumoral NSCLC B cells
#' BCR <- read.csv("https://ncborcherding.github.io/vignettes/b_contigs.csv", 
#' stringsAsFactors = FALSE)
#' combined <- combineBCR(BCR, samples = "Patient1", ID = "Time1")
#' 
#' @param df List of filtered contig annotations from 10x Genomics.
#' @param samples The labels of samples.
#' @param ID The additional sample labeling option.
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @import dplyr
#' @importFrom Biostrings stringDist
#' @export
#' @return List of clonotypes for individual cell barcodes
combineBCR <- function(df, samples = NULL, ID = NULL, removeNA = FALSE, 
                        removeMulti = FALSE) {
    df <- checkList(df)
    out <- NULL
    final <- list()
    checkContigBarcodes(df, samples, ID)
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain %in% c("IGH", "IGK", "IGL"))
        df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        df[[i]] <- df[[i]] %>% group_by(barcode,chain) %>% top_n(n=1,wt=reads)
        df[[i]]$sample <- samples[i]
        df[[i]]$ID <- ID[i]
        df[[i]] <- filteringMulti(df[[i]]) }
    out <- modifyBarcodes(df, samples, ID)
    for (i in seq_along(out)) { 
        data2 <- data.frame(out[[i]])
        data2 <- makeGenes(cellType = "B", data2)
        unique_df <- unique(data2$barcode)
        Con.df <- data.frame(matrix(NA, length(unique_df), 9))
        colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
        Con.df$barcode <- unique_df
        Con.df <- parseBCR(Con.df, unique_df, data2)
        Con.df <- assignCT(cellType = "B", Con.df)
        data3<-Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>%
            mutate(length2 = nchar(cdr3_nt2))
        final[[i]] <- data3 }
    dictionary <- bind_rows(final)
    IGH <- hammingCompare(dictionary, "IGH", "cdr3_nt1", "length1")
    IGLC <- hammingCompare(dictionary, "IGLC", "cdr3_nt2", "length2")
    for(i in seq_along(final)) {
        final[[i]]<-merge(final[[i]],IGH,by.x="cdr3_nt1",by.y="IG",all.x=TRUE)
        final[[i]]<-merge(final[[i]],IGLC,by.x="cdr3_nt2",by.y="IG",all.x=TRUE)
        final[[i]][final[[i]] == "NA_NA" | final[[i]] == "NA_NA_NA_NA"] <- NA 
        num <- ncol(final[[i]])
        final[[i]][,"CTstrict"] <- paste0(final[[i]][,num-1],"_",
        final[[i]][,"vgene1"],"_",final[[i]][,num],"_",final[[i]][,"vgene2"])
        final[[i]]$cellType <- "B"
        final[[i]]$sample <- samples[i]
        final[[i]]$ID <- ID[i]
        final[[i]]<- final[[i]][, c("barcode", "sample", "ID", 
            heavy_lines[c(1,2,3)], light_lines[c(1,2,3)], CT_lines)]}
    names <- NULL
    for (i in seq_along(samples)){c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names
    for (i in seq_along(final)) {
        final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]}
    if (removeNA == TRUE) { final <- removingNA(final) }
    if (removeMulti == TRUE) { final <- removingMulti(final) }
    return(final) }

# Calculates the normalized Hamming Distance between the contig 
# nucleotide sequence.
#' @importFrom Biostrings stringDist
hammingCompare <- function(Con.df, gene, chain, length) {
    `%!in%` = Negate(`%in%`)
    overlap <- NULL
    out <- NULL
    lengths_IGH<-Con.df[duplicated(Con.df[,"length1"]),][,"length1"]
    lengths_IGH <- na.omit(unique(lengths_IGH))
    lengths_IGL<-Con.df[duplicated(Con.df[,"length2"]),][,"length2"]
    lengths_IGL <- na.omit(unique(lengths_IGL))
    specificLength  <- if(gene=="IGH") lengths_IGH else lengths_IGL
    for (i in seq_along(lengths_IGH)) {
        tmp<-na.omit(Con.df[Con.df[,length] == specificLength[i],])
        tmp2 <- as.matrix(stringDist(tmp[,chain], 
                    method = "hamming")/specificLength[i])
        filtered <- which(tmp2 >= 0.85, arr.ind = TRUE)
        if (nrow(filtered) == 0) { next()
        } else if (nrow(filtered) != 0) {
            for (x in seq_along(nrow(filtered))) {
                df <- c(tmp[,chain][filtered[x,1]],tmp[,chain][filtered[x,2]])
                df <- df[order(df)]
                out <- rbind.data.frame(out,df, stringsAsFactors = FALSE)
                out <- unique(out)
                out <- as.data.frame(out, stringsAsFactors = FALSE) }}
        overlap <- rbind.data.frame(overlap,out, stringsAsFactors = FALSE) }
    if (!is.null(overlap)) { colnames(overlap) <- c("Col1", "Col2")
        overlap <- unique(overlap)
        IG <- Con.df[Con.df[,chain] %!in% overlap[,1],]
        IG <- na.omit(unique(IG[IG[,chain] %!in% overlap[,2],][,chain]))
        Hclonotype <- paste0(gene, seq_len(length(IG)))
        IG <- data.frame(IG, Hclonotype)
        unique_over <- data.frame(unique(overlap$Col1), 
                            stringsAsFactors = FALSE)
        unique_over$Hclonotype <- paste0(gene, ":HD", ".", 
                                    seq_len(length(unique_over)))
        colnames(unique_over)[1] <- "barcodes"
        overlap <- merge(overlap, unique_over, by.x="Col1", by.y="barcodes")
        barcodeOverlap <- data.frame(unique(c(overlap[,1], overlap[,2])))
        barcodeOverlap$Hclonotype <- NULL
        for (y in seq_along(nrow(barcodeOverlap))) {
            if (barcodeOverlap[y,1] %in% overlap[,"Col1"]) {
                x <- which(overlap[,"Col1"] == barcodeOverlap[y,1])
            }else if (barcodeOverlap[y,1] %in% overlap[,"Col2"]) {
                x <- which(overlap[,"Col2"] == barcodeOverlap[y,1]) }
            x <- x[1]
            barcodeOverlap[y,2] <- overlap[x,"Hclonotype"] }
        colnames(barcodeOverlap) <- colnames(IG)
        IG <- rbind.data.frame(IG, barcodeOverlap, stringsAsFactors = FALSE)
    } else { IG <- Con.df[,chain]
        IG <- na.omit(unique(IG))
        Hclonotype <- paste0(gene, ".", seq_len(length(IG)))
        IG <- data.frame(IG, Hclonotype) } }
