#' Combining the list of T or B cell Receptors from TRUST4 pipeline
#'
#' This function consolidates a list of TCR/BCR sequencing results to the level of 
#' the individual cell barcodes using the same approach as 
#'  \code{\link{combineTCR}} and \code{\link{combineBCR}}. 
#' Using the samples and ID parameters, the function will add the strings 
#' as prefixes to prevent issues with repeated 
#' barcodes. The resulting new barcodes will need to match the seurat or SCE 
#' object in order to use, \code{\link{combineExpression}}. Several 
#' levels of filtering exist - removeMulti are parameters that 
#' control how  the function  deals with barcodes with multiple chains 
#' recovered. Please \href{https://github.com/liulab-dfci/TRUST4}{read more}
#' and cite the TRUST4 pipeline if using this function.
#' 
#' 
#' @examples
#' \dontrun{ 
#' combineTRUST4(contig_list, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' }
#' 
#' @param df List of Contig outputs from TRUST4
#' @param samples The labels of samples (required).
#' @param ID The additional sample labeling (optional).
#' @param cells The type of cell - T cell-AB or T cell-GD, or B cell
#' @param removeNA This will remove any chain without values.
#' @param threshold If combining B cells - the normalized edit 
#' distance to consider. The higher the number the more similarity 
#' of sequence will be used for clustering.
#' @import dplyr
#' @export
#' @return List of clonotypes for individual cell barcodes


combineTRUST4 <- function(df, samples = NULL, ID = NULL, 
                          cells = c("T-AB", "T-GD", "B"), 
                          removeNA = FALSE, 
                          threshold = 0.85) {
    df <- checkList(df)
    out <- NULL
    final <- NULL
    
    chain1 <- cellT(cells)[[1]]
    chain2 <- cellT(cells)[[2]]
    if (cells == "B") {
      chain1 <- "heavy"
      chain2 <- "light"
    }
    cellType <- cellT(cells)[[3]]
    if (cells != "B") {
        cells <- unname(c("T-AB" = "abT", "T-GD" = "gdT")[cells])
    }
    if(cells != "B") {
        for (i in seq_along(df)) {
            colnames(df[[i]])[1] <- "barcode"
            df[[i]] <- df[[i]][df[[i]]$cell_type == cells,]
            df[[i]][df[[i]] == "*"] <- NA
            
            TCRB <- str_split(df[[i]]$chain1, ",", simplify = TRUE)[,seq_len(6)]
            TCRB[TCRB == "*"] <- "None"
            colnames(TCRB) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
            TCRB <- data.frame(barcode = df[[i]][,1], chain = chain2, TCRB)
            TCRA <- str_split(df[[i]]$chain2, ",", simplify = TRUE)[,seq_len(6)]
            TCRA[TCRA == "*"] <- "None"
            colnames(TCRA) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
            TCRA <- data.frame(barcode = df[[i]][,1], chain = chain1, TCRA)
            data2 <- rbind(TCRA, TCRB)
            data2$sample <- samples[i]
            data2$ID <- ID[i]
            data2[data2 == ""] <- NA
            data2[data2 == "NA"] <- NA
            data2 <- makeGenes(cellType, data2, chain1, chain2)
            
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 7))
            colnames(Con.df) <- c("barcode",tcr1_lines, tcr2_lines)
            Con.df$barcode <- unique_df
            Con.df <- parseTCR(Con.df, unique_df, data2)
            Con.df <- assignCT(cellType, Con.df)
            Con.df$cellType <- cellType
            Con.df[Con.df == "NA_NA" | Con.df == "NA_NA_NA_NA"] <- NA 
            data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))], 
                           Con.df, by = "barcode")
            if (!is.null(sample) & !is.null(ID)) {
              data3<-data3[,c("barcode","sample","ID",tcr1_lines,tcr2_lines,
                              CT_lines)] }
            else if (!is.null(sample) & is.null(ID)) {
              data3<-data3[,c("barcode","sample",tcr1_lines,tcr2_lines,
                              CT_lines)] 
            }
            final[[i]] <- data3
        }
    }
    else {
        for (i in seq_along(df)) {
            colnames(df[[i]])[1] <- "barcode"
            df[[i]] <- df[[i]][df[[i]]$cell_type == cells,]
            df[[i]][df[[i]] == "*"] <- NA
            IGH <- str_split(df[[i]]$chain1, ",", simplify = TRUE)[,seq_len(6)]
            IGH[IGH == "*"] <- "None"
            colnames(IGH ) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
            IGH <- data.frame(barcode = df[[i]][,1], chain = chain1, IGH )
            
            IGL <- str_split(df[[i]]$chain2, ",", simplify = TRUE)[,seq_len(6)]
            IGL[IGL == "*"] <- "None"
            colnames(IGL) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
            IGL <- as.data.frame(IGL)
            IGL<- data.frame(barcode = df[[i]][,1], chain = substr(IGL$v_gene, 1,3), IGL)
            data2 <- rbind(IGH, IGL)
            data2$sample <- samples[i]
            data2$ID <- ID[i]
            data2[data2 == ""] <- NA
            data2[data2 == "NA"] <- NA
            
            data2 <- makeGenes(cellType = "B", data2)
            unique_df <- unique(data2$barcode)
            Con.df <- data.frame(matrix(NA, length(unique_df), 9))
            colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
            Con.df$barcode <- unique_df
            Con.df <- parseBCR(Con.df, unique_df, data2)
            Con.df <- assignCT(cellType = "B", Con.df)
            Con.df[Con.df== "NA"] <- NA
            data3<-Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>%
                mutate(length2 = nchar(cdr3_nt2))
            final[[i]] <- data3 
        }
        dictionary <- bind_rows(final)
        IGH <- lvCompare(dictionary, "IGH", "cdr3_nt1", threshold)
        IGLC <- lvCompare(dictionary, "IGLC", "cdr3_nt2", threshold)
        for(i in seq_along(final)) {
            final[[i]]<-merge(final[[i]],IGH,by.x="cdr3_nt1",by.y="IG",all.x=TRUE)
            final[[i]]<-merge(final[[i]],IGLC,by.x="cdr3_nt2",by.y="IG",all.x=TRUE)
            num <- ncol(final[[i]])
            final[[i]][,"CTstrict"] <- paste0(final[[i]][,num-1],"_",
                                              final[[i]][,"vgene1"],"_",final[[i]][,num],"_",final[[i]][,"vgene2"])
            final[[i]]$cellType <- "B"
            final[[i]]$sample <- samples[i]
            final[[i]]$ID <- ID[i]
            final[[i]][final[[i]] == "NA_NA" | final[[i]] == "NA_NA_NA_NA"] <- NA 
            if (!is.null(sample) & !is.null(ID)) {
              final[[i]]<- final[[i]][, c("barcode", "sample", "ID", 
                                          heavy_lines[c(1,2,3)], light_lines[c(1,2,3)], CT_lines)]
            }
            else if (!is.null(sample) & is.null(ID)) {
              final[[i]]<- final[[i]][, c("barcode", "sample", 
                                          heavy_lines[c(1,2,3)], light_lines[c(1,2,3)], CT_lines)]
            }
        }
    }
    if (!is.null(samples)) {
      final <- modifyBarcodes(final, samples, ID)
    } 
    names <- NULL
    for (i in seq_along(samples)){
      if (!is.null(sample) & !is.null(ID)) {
        c <- paste(samples[i], "_", ID[i], sep="")
      } else if (!is.null(sample) & is.null(ID)) {
        c <- paste(samples[i], sep="")
      }
      names <- c(names, c)
    }
    names(final) <- names
    for (i in seq_along(final)){
        final[[i]]<-final[[i]][!duplicated(final[[i]]$barcode),]
        final[[i]]<-final[[i]][rowSums(is.na(final[[i]])) < 10, ]}
    if (removeNA == TRUE) { final <- removingNA(final)}
    return(final) 
}
