#' Combining the list of T Cell Receptor contigs
#'
#' This function consolidates a list of TCR sequencing results to the level of the indiviudal cell barcodes. Using the samples
#' and ID parameters, the function will add the strings as prefixes to prevent issues with repeated barcodes. The resulting
#' new barcodes will need to match the seurat or SCE object in order to use, @seealso \code{\link{combineExpression}}.
#' Several levels of filtering exist - remove or filterMulti are parameters that control how the function deals with barcodes
#' with multiple chains recovered.
#' @examples
#' combineTCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), rep(c("P", "T"), 3), cells ="T-AB")
#' 
#' @param df List of filtered contig annotations from 10x Genomics.
#' @param samples The labels of samples.
#' @param ID The additional sample labeling option.
#' @param cells The type of T cell - T cell-AB or T cell-GD
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 2 corresponding chains with the
#' highest expression for a single barcode.
#' @import dplyr
#' @export
#' @return List of clonotypes for individual cell barcodes
combineTCR <- function(df,
                           samples = NULL,
                           ID = NULL,
                           cells = c("T-AB", "T-GD"),
                           removeNA = FALSE,
                           removeMulti = FALSE,
                            filterMulti = FALSE) {
    df <- if(is(df)[1] != "list") list(df) else df
    out <- NULL
    final <- NULL
    count <- length(unlist(strsplit(df[[1]]$barcode[1], "[-]")))
    count2 <- length(unlist(strsplit(df[[1]]$barcode[1], "[_]")))
    if (count > 2 | count2 > 2) {
        stop("Seems to be an error in the naming of the contigs, ensure the barcodes are labeled like, AAACGGGAGATGGCGT-1 or AAACGGGAGATGGCGT, use stripBarcode to get the basic format", call. = FALSE)
    } else if (length(df) != length(samples) | length(df) != length(ID)) {
        stop("Make sure the sample and ID labels match the length of the list of data frames (df).", call. = FALSE)
    } else {
            if (cells == "T-AB") {
                chain1 <- "TRA"
                chain2 <- "TRB"
                cellType <- "T-AB"
            }
            else if (cells == "T-GD") {
                chain1 <- "TRG"
                chain2 <- "TRD"
                cellType <- "T-GD"
            }
            for (i in seq_along(df)) {
                df[[i]] <- subset(df[[i]], chain != "Multi")
                df[[i]] <- subset(df[[i]], chain == chain1 | chain == chain2)
                df[[i]] <- subset(df[[i]], productive == TRUE | productive == "TRUE" | productive == "True")
                df[[i]]$sample <- samples[i]
                df[[i]]$ID <- ID[i]
                if (filterMulti == TRUE) {
                    barcodes <- as.character(unique(table$Var1))
                    multichain <- NULL
                    for (j in seq_along(barcodes)) {
                        chain <- df[[i]][df[[i]]$barcode == barcodes[j],] %>% group_by(barcode) %>% top_n(n = 2, wt = reads)
                        multichain <- rbind(multichain, chain)
                    }
                    `%!in%` = Negate(`%in%`)
                    df[[i]] <- subset(df[[i]], barcode %!in% barcodes)
                    df[[i]] <- rbind(df[[i]], multichain)
                }
                if (nrow(df[[i]]) == 0) {
                    stop("Check some hypotenuses, Captain. There are 0 contigs after filtering for celltype.", call. = FALSE)
                }
            }

        for (x in seq_along(df)) {
            data <- df[[x]]
            data$barcode <- paste(samples[x], "_", ID[x], "_", data$barcode, sep="")
            out[[x]] <- data
        }

        for (i in seq_along(out)) {

            data2 <- out[[i]]
            data2 <- data2 %>%
                mutate(TCR1 = ifelse(chain == chain1, paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
                mutate(TCR2 = ifelse(chain == chain2, paste(with(data2, interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 7))

            colnames(Con.df) <- c("barcode","TCR1", "cdr3_aa1", "cdr3_nt1", "TCR2", "cdr3_aa2", "cdr3_nt2")
            Con.df$barcode <- unique_df
            y <- NULL
            for (y in seq_along(unique_df)){
                barcode.i <- Con.df$barcode[y]
                location.i <- which(barcode.i == data2$barcode)
                if (length(location.i) == 2){
                    if (is.na(data2[location.i[1],c("TCR1")])) {
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")]
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[2],c("TCR1", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")]
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")]
                    }
                } else if (length(location.i) == 3) { # if there are more than 2 chains for each cell
          if (is.na(data2[location.i[1],c("TCR1")])) { # if the 1st location is occupied by TRB
            Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")] # add TRB info to the TCR2 related columns
            if (is.na(data2[location.i[2],c("TCR1")])) { # if the 2nd location is occupied by TRB
              TRdf <- paste(Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")],data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRBs from 1st and 2nd locations
              Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add the combination of 2 TRB chains to the TCR2 slot
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")] # add TRA from the 3rd location to the TCR1 related columns
            } else { # if the 2nd location is occupied by TRA
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")] # add TRA info to the TCR1 related columns
              if (is.na(data2[location.i[3],c("TCR1")])) { # if the 3rd location is occupied by TRB
                TRdf <- paste(Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")],data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRB from 1st and 3rd locations
                Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add combined TRB relatedn info to the TCR2 related columns
              } else { # if the 3rd location is occupied by TRA
                TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRAs from 2nd and 3rd locations
                Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # Add combined TRA to TCR1 related columns
              }
            }
          } else { # if 1st location is occupied by TRA
            Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")] # add 1st TRA to TCR1 column
            if (is.na(data2[location.i[2],c("TCR1")])) { # if the 2nd location is occupied by TRB
              if (is.na(data2[location.i[3],c("TCR1")])) { # if the 3rd location is occupied by TRB
                TRdf <- paste(data2[location.i[2],c("TCR2", "cdr3", "cdr3_nt")],data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")],sep=";") # paste TRB from 2nd and 3rd locations
                Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- TRdf # add TRB combo to TCR2 related  columns
              } else { # if TRA is on 3rd location
                TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[3],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRA from 1st and 3rd locations
                Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # add TRA combo to TCR1 related columns
              }
            } else { # if TRA is on 2nd location
              TRdf <- paste(Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")],data2[location.i[2],c("TCR1", "cdr3", "cdr3_nt")],sep=";") # paste TRA from 1st and 2nd locations
              Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- TRdf # add the combination of 2 TRA chains to the TCR1 slot
              Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[3],c("TCR2", "cdr3", "cdr3_nt")] # add TRB to the TCR2 related columns
            }
          }
        } else if (length(location.i) == 1) {
                    chain.i <- data2$chain[location.i]
                    if (chain.i == "TRA"){
                        Con.df[y,c("TCR1", "cdr3_aa1", "cdr3_nt1")] <- data2[location.i[1],c("TCR1", "cdr3", "cdr3_nt")]
                    } else {
                        Con.df[y,c("TCR2", "cdr3_aa2", "cdr3_nt2")] <- data2[location.i[1],c("TCR2", "cdr3", "cdr3_nt")]
                    }
                }
            }
            Con.df$CTgene <- paste(Con.df$TCR1, Con.df$TCR2, sep="_")
            Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
            Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
            Con.df$CTstrict <- paste(Con.df$TCR1, Con.df$cdr3_nt1, Con.df$TCR2, Con.df$cdr3_nt2, sep="_")
            Con.df$cellType <- cells #autodetect chains in new functions
            Con.df[Con.df == "NA_NA"] <- NA #remove the na when gene, aa, or nt is called later
            Con.df[Con.df == "NA_NA_NA_NA"] <- NA #remove the na when nt+gene is called later
            data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))], Con.df, by = "barcode")
            data3 <- data3[, c("barcode", "sample", "ID", "TCR1", "cdr3_aa1", "cdr3_nt1", "TCR2", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", "CTstrict", "cellType")]
            final[[i]] <- data3
        }

    names <- NULL
    for (i in seq_along(samples)) {
        c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names

    for (i in seq_along(final)) {
    final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]
    }
    if (removeNA == TRUE) {
        for(i in seq_along(final)) {
            final[[i]] <- na.omit(final[[i]])
        }
    }
    if (removeMulti == TRUE) {
        for(i in seq_along(final)) {
            final[[i]] <- filter(final[[i]], !grepl(";",CTnt))
        }
    }
    }
    return(final)
}

#' Combining the list of B Cell Receptor contigs
#'
#' This function consolidates a list of BCR sequencing results to the level of the indiviudal cell barcodes. Using the samples
#' and ID parameters, the function will add the strings as prefixes to prevent issues with repeated barcodes. The resulting
#' new barcodes will need to match the seurat or SCE object in order to use, @seealso \code{\link{combineExpression}}.
#' Unlike combineTCR(), combineBCR produces a column CTstrict of an index of nucleotide sequence and the corresponding v-gene.
#' This index automatically caluclates the Hammings distance between sequences of the same length and will index sequences with
#' > 0.85 normalized Hammings distance with the same ID. If nucleotide sequences meet the threshold, ":HD" will be added to
#' the CTstrict column string.
#'
#' @examples
#' \donttest{
#' combineBCR(contig_list, rep(c("PX", "PY", "PZ"), each=2), rep(c("P", "T"), 3))
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
combineBCR <- function(df,
                       samples = NULL,
                       ID = NULL,
                       removeNA = FALSE,
                       removeMulti = FALSE) {
    df <- if(is(df)[1] != "list") list(df) else df
    out <- NULL
    final <- list()
    count <- length(unlist(strsplit(df[[1]]$barcode[1], "[-]")))
    count2 <- length(unlist(strsplit(df[[1]]$barcode[1], "[_]")))
    if (count > 2 | count2 > 2) {
        stop("Seems to be an error in the naming of the contigs, ensure the barcodes are labeled like, AAACGGGAGATGGCGT-1 or AAACGGGAGATGGCGT, use stripBarcode to get the basic format", call. = FALSE)
    } else if (length(df) != length(samples) | length(df) != length(ID)) {
        stop("Make sure the sample and ID labels match the length of the list of data frames (df).", call. = FALSE)
    }
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain == "IGH" | chain == "IGK" | chain == "IGL")
        df[[i]] <- subset(df[[i]], productive == TRUE | productive == "TRUE" | productive == "True")
        df[[i]] <- df[[i]] %>%
            group_by(barcode, chain) %>%
            top_n(n = 1, wt = reads)
        df[[i]]$sample <- samples[i]
        df[[i]]$ID <- ID[i]
        table <- subset(as.data.frame(table(df[[i]]$barcode)), Freq > 2)
        barcodes <- as.character(unique(table$Var1))
        multichain <- NULL
        for (j in seq_along(barcodes)) {
            chain <- df[[i]][df[[i]]$barcode == barcodes[j],] %>% group_by(barcode) %>% top_n(n = 2, wt = reads)
            multichain <- rbind(multichain, chain)
        }
        `%!in%` = Negate(`%in%`)
        df[[i]] <- subset(df[[i]], barcode %!in% barcodes)
        df[[i]] <- rbind(df[[i]], multichain)
    }
    for (x in seq_along(df)) {
        data <- df[[x]]
        data$barcode <- paste(samples[x], "_", ID[x], "_", data$barcode, sep="")
        out[[x]] <- data
    }

    for (i in seq_along(out)) {
        data2 <- data.frame(out[[i]], stringsAsFactors = FALSE)
        data2 <- data2 %>%
            mutate(IGKct = ifelse(chain == "IGK", paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGLct = ifelse(chain == "IGL", paste(with(data2, interaction(v_gene,  j_gene, c_gene))), NA)) %>%
            mutate(IGHct = ifelse(chain == "IGH", paste(with(data2, interaction(v_gene,  j_gene, d_gene, c_gene))), NA))
        unique_df <- unique(data2$barcode)
        Con.df <- data.frame(matrix(NA, length(unique_df), 9 ))
        colnames(Con.df) <- c("barcode","IGH", "cdr3_aa1", "cdr3_nt1", "vgene1", "IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")
        Con.df$barcode <- unique_df
        y <- NULL
        for (y in seq_along(unique_df)){
            barcode.i <- Con.df$barcode[y]
            location.i <- which(barcode.i == data2$barcode)
            if (length(location.i) == 2){
                if (is.na(data2[location.i[1],c("IGHct")])) {
                    Con.df[y,c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")] <- data2[location.i[2],c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                    if (is.na(data2[location.i[1],c("IGKct")])) {
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[1],c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
                    }
                    else {
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[1],c("IGKct", "cdr3", "cdr3_nt", "v_gene")]
                    }
                } else {
                    Con.df[y,c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")] <- data2[location.i[1],c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                    if (is.na(data2[location.i[1],c("IGKct")])) {
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[2],c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
                    }
                    else {
                        Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[1],c("IGKct", "cdr3", "cdr3_nt", "v_gene")]
                    }
                }
            } else if (length(location.i) == 1) {
                chain.i <- data2$chain[location.i]
                if (chain.i == "IGH"){
                    Con.df[y,c("IGH", "cdr3_aa1", "cdr3_nt1", "vgene1")] <- data2[location.i[1],c("IGHct", "cdr3", "cdr3_nt", "v_gene")]
                } else if (chain.i == "IGL") {
                    Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[2],c("IGLct", "cdr3", "cdr3_nt", "v_gene")]
                }
                else {
                    Con.df[y,c("IGLC", "cdr3_aa2", "cdr3_nt2", "vgene2")] <- data2[location.i[1],c("IGKct", "cdr3", "cdr3_nt", "v_gene")]
                }
            }
        }

    Con.df$CTgene <- paste(Con.df$IGH, Con.df$IGLC, sep="_")
    Con.df$CTnt <- paste(Con.df$cdr3_nt1, Con.df$cdr3_nt2, sep="_")
    Con.df$CTaa <- paste(Con.df$cdr3_aa1, Con.df$cdr3_aa2, sep="_")
    data3 <- Con.df %>%
        mutate(length1 = nchar(cdr3_nt1)) %>%
        mutate(length2 = nchar(cdr3_nt2))
    final[[i]] <- data3
    }
    #Need to compute normalized hamming distance for across all nucleotide sequences in both IGH and the light chain
    dictionary <- bind_rows(final)
    IGH <- hammingCompare(dictionary, "IGH", "cdr3_nt1", "length1")
    IGLC <- hammingCompare(dictionary, "IGLC", "cdr3_nt2", "length2")

    for(i in seq_along(final)) {
        final[[i]] <- merge(final[[i]], IGH, by.x = "cdr3_nt1", by.y = "IG", all.x = TRUE)
        final[[i]] <- merge(final[[i]], IGLC, by.x = "cdr3_nt2", by.y = "IG", all.x = TRUE)
        final[[i]][,"CTstrict"] <- paste0(final[[i]][,ncol(final[[i]])-1], "_", final[[i]][,"vgene1"], "_", final[[i]][,ncol(final[[i]])], "_", final[[i]][,"vgene2"])
        final[[i]]$cellType <- "B"
        final[[i]]$sample <- samples[i]
        final[[i]]$ID <- ID[i]
        final[[i]]<- final[[i]][, c("barcode", "sample", "ID", "IGH", "cdr3_aa1", "cdr3_nt1", "IGLC", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", "CTstrict", "cellType")]
    }

    names <- NULL
    for (i in seq_along(samples)) {
        c <- paste(samples[i], "_", ID[i], sep="")
        names <- c(names, c)}
    names(final) <- names

    for (i in seq_along(final)) {
        final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]
    }
    if (removeNA == TRUE) {
        for(i in seq_along(final)) {
            final[[i]] <- na.omit(final[[i]])
        }
    }
    if (removeMulti == TRUE) {
        for(i in seq_along(final)) {
            final[[i]] <- filter(final[[i]], !grepl(";",CTnt))
        }
    }
    return(final)
}

#' Calculates the normalized Hamming Distance between the contig nucleotide sequence.
#'
#' This feature caluclates the normalized Hammings Distance, or the Hammings Distance divided by length of
#' sequence in order to index similar and divergent sequences involved in combineBCR(). It is not designed
#' as an indepenent function. The threshold for similar sequences is set to > 0.85 normalized Hammings
#' Distance, which if met will index the similar sequence into a single sequence and add ":HD" to the index.
#' @keywords internal
#' @param gene The IGH or IG light chains (IGLC)
#' @param chain The column header with the nucletoide sequence
#' @param length The column header with the specific length
#' @importFrom Biostrings stringDist
#' @return #Normalized Hammings Distance comparisons to be used in combineBCR()
hammingCompare <- function(Con.df, gene, chain, length) {
    `%!in%` = Negate(`%in%`)
    overlap <- NULL
    out <- NULL
    lengths_IGH <- Con.df[duplicated(Con.df[,"length1"]),]
    lengths_IGH <- na.omit(unique(lengths_IGH[,"length1"]))
    lengths_IGL <- Con.df[duplicated(Con.df[,"length2"]),]
    lengths_IGL <- na.omit(unique(lengths_IGL[,"length2"]))
    if (gene == "IGH") {
        specificLength <- lengths_IGH
    } else if (gene == "IGLC") {
        specificLength <- lengths_IGL
    }
    for (i in seq_along(lengths_IGH)) {
        tmp <- na.omit(Con.df[Con.df[,length] == specificLength[i],])
        tmp2 <- as.matrix(stringDist(tmp[,chain], method = "hamming")/specificLength[i])
        filtered <- which(tmp2 >= 0.85, arr.ind = TRUE)
        vec <- nrow(filtered)
        if (vec == 0) {
            next()
        } else if (vec != 0) {
            for (x in seq_along(vec)) {
                df <- c(tmp[,chain][filtered[x,1]], tmp[,chain][filtered[x,2]])
                df <- df[order(df)]
                out <- rbind.data.frame(out,df, stringsAsFactors = FALSE)
                out <- unique(out)
                out <- as.data.frame(out, stringsAsFactors = FALSE)
            }
        }
        overlap <- rbind.data.frame(overlap,out, stringsAsFactors = FALSE)
    }
    if (!is.null(overlap)) {
        colnames(overlap) <- c("Col1", "Col2")
        overlap <- unique(overlap)

        IG <- Con.df[Con.df[,chain] %!in% overlap[,1],]
        IG <- IG[IG[,chain] %!in% overlap[,2],]
        IG <- na.omit(unique(IG[,chain]))
        Hclonotype <- paste0(gene, seq_len(IG))
        IG <- data.frame(IG, Hclonotype)
        unique_over <- data.frame(unique(overlap$Col1), stringsAsFactors = FALSE)
        unique_over$Hclonotype <- paste0(gene, ":HD", ".", seq_len(unique_over))
        colnames(unique_over)[1] <- "barcodes"
        overlap <- merge(overlap, unique_over, by.x="Col1", by.y="barcodes")
        barcodeOverlap <- unique(c(overlap[,1], overlap[,2]))
        barcodeOverlap <- data.frame(barcodeOverlap, stringsAsFactors = FALSE)
        barcodeOverlap$Hclonotype <- NULL
        vec <- nrow(barcodeOverlap)
        for (y in seq_along(vec)) {
            if (barcodeOverlap[y,1] %in% overlap[,"Col1"]) {
                x <- which(overlap[,"Col1"] == barcodeOverlap[y,1])
            }else if (barcodeOverlap[y,1] %in% overlap[,"Col2"]) {
                x <- which(overlap[,"Col2"] == barcodeOverlap[y,1])
            }
            x <- x[1]
            barcodeOverlap[y,2] <- overlap[x,"Hclonotype"]
        }
        colnames(barcodeOverlap) <- colnames(IG)
        IG <- rbind.data.frame(IG, barcodeOverlap, stringsAsFactors = FALSE)
    }
    else {
        IG <- Con.df[,chain]
        IG <- na.omit(unique(IG))
        Hclonotype <- paste0(gene, ".", seq_len(IG))
        IG <- data.frame(IG, Hclonotype)
    }
}
