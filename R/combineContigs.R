#' 
#' Combining the list of T Cell Receptor contigs
#'
#' This function consolidates a list of TCR sequencing results to the level of 
#' the individual cell barcodes. Using the samples and ID parameters, the 
#' function will add the strings as prefixes to prevent issues with repeated 
#' barcodes. The resulting new barcodes will need to match the seurat or SCE 
#' object in order to use, @seealso \code{\link{combineExpression}}. Several 
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
                        cells = c("T-AB", "T-GD"), 
                        removeNA = FALSE, removeMulti = FALSE, 
                        filterMulti = FALSE) {
    df <- checkList(df)
    tcr1_lines <- c("TCR1", "cdr3_aa1", "cdr3_nt1")
    tcr2_lines <- c("TCR2", "cdr3_aa2", "cdr3_nt2")
    data1_lines <- c("TCR1", "cdr3", "cdr3_nt")
    data2_lines <- c("TCR2", "cdr3", "cdr3_nt")
    out <- NULL
    final <- NULL
    checkContigBarcodes(df, samples, ID)
        if (cells == "T-AB") { chain1 <- "TRA" 
            chain2 <- "TRB"
            cellType <- "T-AB"
        } else if (cells == "T-GD") { chain1 <- "TRG" 
            chain2 <- "TRD"
            cellType <- "T-GD" }
        for (i in seq_along(df)) {
            df[[i]] <- subset(df[[i]], chain != "Multi")
            df[[i]] <- subset(df[[i]], chain == chain1 | chain == chain2)
            df[[i]] <- subset(df[[i]], productive == TRUE | 
                                productive == "TRUE" | productive == "True")
            df[[i]]$sample <- samples[i]
            df[[i]]$ID <- ID[i]
            if (filterMulti == TRUE) {
                barcodes <- as.character(unique(table$Var1))
                multichain <- NULL
                for (j in seq_along(barcodes)) {
                    chain <- df[[i]][df[[i]]$barcode == barcodes[j],] %>% 
                            group_by(barcode) %>% top_n(n = 2, wt = reads)
                    multichain <- rbind(multichain, chain) }
                `%!in%` = Negate(`%in%`)
                df[[i]] <- subset(df[[i]], barcode %!in% barcodes)
                df[[i]] <- rbind(df[[i]], multichain) }
            if (nrow(df[[i]]) == 0) {
                stop("Check some hypotenuses, Captain. There are 0 contigs 
                    after filtering for celltype.", call. = FALSE) }}
        out <- modifyBarcodes(df, samples, ID)
        for (i in seq_along(out)) { data2 <- out[[i]]
            data2 <- data2 %>% mutate(TCR1 = ifelse(chain == chain1, 
                    paste(with(data2, 
                    interaction(v_gene,  j_gene, c_gene))), NA)) %>%
                mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, 
                    interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 7))
            colnames(Con.df) <- c("barcode",tcr1_lines, tcr2_lines)
            Con.df$barcode <- unique_df
            for (y in seq_along(unique_df)){
                barcode.i <- Con.df$barcode[y]
                location.i <- which(barcode.i == data2$barcode)
                if (length(location.i) == 2){
                    if (is.na(data2[location.i[1],c("TCR1")])) {
                        Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines]
                        Con.df[y,tcr1_lines]<-data2[location.i[2],data1_lines]
                    } else {
                        Con.df[y,tcr1_lines]<-data2[location.i[1],data1_lines]
                        Con.df[y,tcr2_lines]<-data2[location.i[2],data2_lines] }
                } else if (length(location.i) == 3) { 
                    if (is.na(data2[location.i[1],c("TCR1")])) { 
                        Con.df[y,tcr2_lines]<-data2[location.i[1],data2_lines] 
                        if (is.na(data2[location.i[2],c("TCR1")])) { 
                            TRdf <- paste(Con.df[y, tcr2_lines],
                                    data2[location.i[2], data2_lines],sep=";") 
                            Con.df[y,tcr2_lines] <- TRdf 
                            Con.df[y,tcr1_lines] <- 
                                data2[location.i[3],data1_lines] 
                        } else { # if the 2nd location is occupied by TRA
                            Con.df[y,tcr1_lines] <- 
                                data2[location.i[1],data1_lines] 
                        if (is.na(data2[location.i[3],c("TCR1")])) { 
                            TRdf <- paste(Con.df[y, tcr2_lines],
                                    data2[location.i[3], data2_lines],sep=";") 
                            Con.df[y,tcr2_lines] <- TRdf 
                        } else { # if the 3rd location is occupied by TRA
                            TRdf <- paste(Con.df[y, tcr1_lines],
                                    data2[location.i[3], data1_lines],sep=";") 
                            Con.df[y,tcr1_lines] <- TRdf }}
                } else { # if 1st location is occupied by TRA
                    Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines] 
                    if (is.na(data2[location.i[2],c("TCR1")])) { 
                        if (is.na(data2[location.i[3],c("TCR1")])) { 
                            TRdf <- paste(data2[location.i[2], data2_lines],
                                    data2[location.i[3], data2_lines],sep=";") 
                            Con.df[y,tcr2_lines] <- TRdf 
                        } else { # if TRA is on 3rd location
                            TRdf <- paste(Con.df[y, tcr1_lines],
                                    data2[location.i[3],data1_lines],sep=";") 
                            Con.df[y,tcr1_lines] <- TRdf }
                    } else { # if TRA is on 2nd location
                        TRdf <- paste(Con.df[y, tcr1_lines],
                                data2[location.i[2], data1_lines],sep=";") 
                        Con.df[y,tcr1_lines] <- TRdf 
                        Con.df[y,tcr2_lines] <- 
                            data2[location.i[3],data2_lines] } }
        } else if (length(location.i) == 1) {
            chain.i <- data2$chain[location.i]
            if (chain.i == "TRA"){
                Con.df[y,tcr1_lines] <- data2[location.i[1],data1_lines]
            } else {Con.df[y,tcr2_lines] <- data2[location.i[1],data2_lines]}}}
        Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
        Con.df$CTstrict <- paste(Con.df$TCR1, Con.df$cdr3_nt1, 
                            Con.df$TCR2, Con.df$cdr3_nt2, sep="_")
        Con.df$cellType <- cells 
        Con.df[Con.df == "NA_NA" | Con.df == "NA_NA_NA_NA"] <- NA 
        data3 <- merge(data2[,-which(names(data2) %in% 
                                c("TCR1","TCR2"))], Con.df, by = "barcode")
        data3 <- data3[, c("barcode", "sample", "ID", "TCR1", 
                            tcr1_lines, tcr2_lines, "CTgene", "CTnt", "CTaa", 
                            "CTstrict", "cellType")]
        final[[i]] <- data3 }
    names <- NULL
    for (i in seq_along(samples)) { c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names
    for (i in seq_along(final)) {
    final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),] }
    if (removeNA == TRUE) { final <- removingNA(final)}
    if (removeMulti == TRUE) { final <- removingMulti(final) }
    return(final)
}

#' Combining the list of B Cell Receptor contigs
#'
#' This function consolidates a list of BCR sequencing results to the level 
#' of the individual cell barcodes. Using the samples and ID parameters, 
#' the function will add the strings as prefixes to prevent issues with 
#' repeated barcodes. The resulting new barcodes will need to match the 
#' seurat or SCE object in order to use, 
#' @seealso \code{\link{combineExpression}}. Unlike combineTCR(), 
#' combineBCR produces a column CTstrict of an index of nucleotide sequence 
#' and the corresponding v-gene. This index automatically caluclates 
#' the Hammings distance between sequences of the same length and will 
#' index sequences with > 0.85 normalized Hammings distance with the same 
#' ID. If nucleotide sequences meet the threshold, ":HD" will be added to 
#' the CTstrict column string.
#'
#' @examples
#' \donttest{
#' combineBCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3))
#' }
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
    heavy_lines <- c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")
    light_lines <- c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")
    out <- NULL
    final <- list()
    checkContigBarcodes(df, samples, ID)
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain == "IGH" | chain == "IGK" 
                    | chain == "IGL")
        df[[i]] <- subset(df[[i]], productive == TRUE 
                    | productive == "TRUE" | productive == "True")
        df[[i]] <- df[[i]] %>% group_by(barcode, chain) %>%
            top_n(n = 1, wt = reads)
        df[[i]]$sample <- samples[i]
        df[[i]]$ID <- ID[i]
        table <- subset(as.data.frame(table(df[[i]]$barcode)), Freq > 2)
        barcodes <- as.character(unique(table$Var1))
        multichain <- NULL
        for (j in seq_along(barcodes)) {
            chain <- df[[i]][df[[i]]$barcode == barcodes[j],] %>% 
                        group_by(barcode) %>% 
                        top_n(n = 2, wt = reads)
            multichain <- rbind(multichain, chain) }
        `%!in%` = Negate(`%in%`)
        df[[i]] <- subset(df[[i]], barcode %!in% barcodes)
        df[[i]] <- rbind(df[[i]], multichain) }
    out <- modifyBarcodes(df, samples, ID)
    for (i in seq_along(out)) {
        data2 <- data.frame(out[[i]], stringsAsFactors = FALSE)
        data2 <- data2 %>%
            mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, 
                interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, 
                interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, 
                interaction(v_gene, j_gene, d_gene, c_gene))), NA))
        unique_df <- unique(data2$barcode)
        Con.df <- data.frame(matrix(NA, length(unique_df), 9 ))
        colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
        Con.df$barcode <- unique_df
        y <- NULL
        for (y in seq_along(unique_df)){
            barcode.i <- Con.df$barcode[y]
            location.i <- which(barcode.i == data2$barcode)
            if (length(location.i) == 2){
                if (is.na(data2[location.i[1],c("IGHct")])) {
                    Con.df[y,heavy_lines] <- data2[location.i[2],
                        c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                    if (is.na(data2[location.i[1],c("IGKct")])) {
                        Con.df[y,light_lines]<- data2[location.i[1],
                        c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
                    } else {
                        Con.df[y,light_lines]<- data2[location.i[1],
                        c("IGKct", "cdr3", "cdr3_nt", "v_gene")]}
                } else { Con.df[y,heavy_lines] <- 
                        data2[location.i[1],
                        c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                    if (is.na(data2[location.i[1],c("IGKct")])) {
                        Con.df[y,light_lines]<- data2[location.i[2],
                        c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
                    } else {
                        Con.df[y,light_lines]<- data2[location.i[1],
                        c("IGKct", "cdr3", "cdr3_nt", "v_gene")]}}
            } else if (length(location.i) == 1) {
                chain.i <- data2$chain[location.i]
                if (chain.i == "IGH"){
                    Con.df[y,heavy_lines]<-data2[location.i[1],
                    c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                } else if (chain.i == "IGL") {
                    Con.df[y,light_lines]<- data2[location.i[2],
                    c("IGLct", "cdr3", "cdr3_nt", "v_gene")]}
                else {
                    Con.df[y,light_lines]<- data2[location.i[1],
                    c("IGKct", "cdr3", "cdr3_nt", "v_gene")] }}}
        Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
        Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
        Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
        data3 <- Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>%
            mutate(length2 = nchar(cdr3_nt2))
        final[[i]] <- data3 }
    dictionary <- bind_rows(final)
    IGH <- hammingCompare(dictionary, "IGH", "cdr3_nt1", "length1")
    IGLC <- hammingCompare(dictionary, "IGLC", "cdr3_nt2", "length2")
    for(i in seq_along(final)) {
        final[[i]] <- merge(final[[i]], IGH, 
                            by.x = "cdr3_nt1", by.y = "IG", all.x = TRUE)
        final[[i]] <- merge(final[[i]], IGLC, 
                            by.x = "cdr3_nt2", by.y = "IG", all.x = TRUE)
        final[[i]][,"CTstrict"] <- paste0(final[[i]][,ncol(final[[i]])-1], 
                                        "_", final[[i]][,"vgene1"], "_", 
                                        final[[i]][,ncol(final[[i]])], "_", 
                                        final[[i]][,"vgene2"])
        final[[i]]$cellType <- "B"
        final[[i]]$sample <- samples[i]
        final[[i]]$ID <- ID[i]
        final[[i]]<- final[[i]][, 
            c("barcode", "sample", "ID", "IGH", "cdr3_aa1", "cdr3_nt1", 
            "IGLC", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", 
            "CTstrict", "cellType")]}
    names <- NULL
    for (i in seq_along(samples)) {
        c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names
    for (i in seq_along(final)) {
        final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]}
    if (removeNA == TRUE) { final <- removingNA(final) }
    if (removeMulti == TRUE) { final <- removingMulti(final) }
    return(final)
}

# Calculates the normalized Hamming Distance between the contig 
# nucleotide sequence.
#' @importFrom Biostrings stringDist
hammingCompare <- function(Con.df, gene, chain, length) {
    `%!in%` = Negate(`%in%`)
    overlap <- NULL
    out <- NULL
    lengths_IGH <- na.omit(unique(Con.df[duplicated(Con.df[,"length1"]),]))
    lengths_IGL <- na.omit(unique(Con.df[duplicated(Con.df[,"length2"]),]))
    if (gene == "IGH") { specificLength <- lengths_IGH
    } else if (gene == "IGLC") { specificLength <- lengths_IGL }
    for (i in seq_along(lengths_IGH)) {
        tmp <- na.omit(Con.df[Con.df[,length] == specificLength[i],])
        tmp2 <- as.matrix(stringDist(tmp[,chain], 
                    method = "hamming")/specificLength[i])
        filtered <- which(tmp2 >= 0.85, arr.ind = TRUE)
        if (nrow(filtered) == 0) { next()
        } else if (nrow(filtered) != 0) {
            for (x in seq_along(nrow(filtered))) {
                df <- c(tmp[,chain][filtered[x,1]], 
                        tmp[,chain][filtered[x,2]])
                df <- df[order(df)]
                out <- rbind.data.frame(out,df, stringsAsFactors = FALSE)
                out <- unique(out)
                out <- as.data.frame(out, stringsAsFactors = FALSE) }}
        overlap <- rbind.data.frame(overlap,out, stringsAsFactors = FALSE) }
    if (!is.null(overlap)) { colnames(overlap) <- c("Col1", "Col2")
        overlap <- unique(overlap)
        IG <- Con.df[Con.df[,chain] %!in% overlap[,1],]
        IG <- na.omit(unique(IG[IG[,chain] %!in% overlap[,2],][,chain]))
        Hclonotype <- paste0(gene, seq_len(IG))
        IG <- data.frame(IG, Hclonotype)
        unique_over <- data.frame(unique(overlap$Col1), 
                            stringsAsFactors = FALSE)
        unique_over$Hclonotype <- paste0(gene, ":HD", ".", 
                                    seq_len(unique_over))
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
        Hclonotype <- paste0(gene, ".", seq_len(IG))
        IG <- data.frame(IG, Hclonotype) } }

# data('v_gene','j_gene', 'c_gene', 'd_gene')
utils::globalVariables(c("v_gene", "j_gene", "c_gene", "d_gene"))
