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
CT_lines <- c("CTgene", "CTnt", "CTaa", "CTstrict")

utils::globalVariables(c(
    "heavy_lines", "light_lines", "l_lines", "k_lines", "h_lines", "tcr1_lines",
    "tcr2_lines", "data1_lines", "data2_lines", "CT_lines"
))

#' @title Combining the list of T cell receptor contigs into clones
#'
#' @description This function consolidates a list of TCR sequencing results to
#' the level of  the individual cell barcodes. Using the \strong{samples} and 
#' \strong{ID} parameters, the function will add the strings as prefixes to 
#' prevent issues with repeated  barcodes. The resulting new barcodes will 
#' need to match the Seurat or SCE object in order to use, 
#' \code{\link{combineExpression}}. Several levels of filtering exist - 
#' \emph{removeNA}, \emph{removeMulti}, or \emph{filterMulti} are parameters 
#' that control how the function deals with barcodes with multiple chains 
#' recovered.
#' 
#' @examples
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' @param input.data List of filtered contig annotations or 
#' outputs from \code{\link{loadContigs}}.
#' @param samples The labels of samples (recommended).
#' @param ID The additional sample labeling (optional).
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 2 
#' corresponding chains with the highest expression for a single barcode. 
#' @param filterNonproductive This option will allow for the removal of 
#' nonproductive chains if the variable exists in the contig data. Default
#' is set to TRUE to remove nonproductive contigs.
#' 
#' @import dplyr
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return List of clones for individual cell barcodes
#' 
combineTCR <- function(input.data, 
                       samples = NULL, 
                       ID = NULL, 
                       removeNA = FALSE, 
                       removeMulti = FALSE, 
                       filterMulti = FALSE,
                       filterNonproductive = TRUE) {
    input.data <- .checkList(input.data)
    input.data <- .checkContigs(input.data)
    out <- NULL
    final <- NULL
    for (i in seq_along(input.data)) {
        if(c("chain") %in% colnames(input.data[[i]])) {
          input.data[[i]] <- subset(input.data[[i]], chain != "Multi")
        }
        if(c("productive") %in% colnames(input.data[[i]]) & filterNonproductive) {
          input.data[[i]] <- subset(input.data[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        }
        input.data[[i]]$sample <- samples[i]
        input.data[[i]]$ID <- ID[i]
        if (filterMulti) { 
          input.data[[i]] <- .filteringMulti(input.data[[i]]) 
        }
    }
    #Prevents error caused by list containing elements with 0 rows
    blank.rows <- which(unlist(lapply(input.data, nrow)) == 0)
    if(length(blank.rows) > 0) {
      input.data <- input.data[-blank.rows]
      if(!is.null(samples)) {
        samples <- samples[-blank.rows]
      }
      if(!is.null(ID)) {
        ID <- ID[-blank.rows]
      }
    }
    if (!is.null(samples)) {
      out <- .modifyBarcodes(input.data, samples, ID)
    } else {
      out <- input.data
    }
    for (i in seq_along(out)) {
        data2 <- .makeGenes(cellType = "T", out[[i]])
        Con.df <- .constructConDfAndParseTCR(data2)
        Con.df <- .assignCT(cellType = "T", Con.df)
        Con.df[Con.df == "NA_NA" | Con.df == "NA;NA_NA;NA"] <- NA 
        data3 <- merge(data2[,-which(names(data2) %in% c("TCR1","TCR2"))], 
            Con.df, by = "barcode")
        if (!is.null(samples) && !is.null(ID)) {
            data3 <- data3[, c("barcode", "sample", "ID", tcr1_lines, tcr2_lines,
                CT_lines)] }
        else if (!is.null(samples) & is.null(ID)) {
          data3<-data3[,c("barcode","sample",tcr1_lines,tcr2_lines,
                          CT_lines)] 
        }
        final[[i]] <- data3 
    }
    name_vector <- character(length(samples))
    for (i in seq_along(samples)) { 
        if (!is.null(samples) && !is.null(ID)) {
            curr <- paste(samples[i], "_", ID[i], sep="")
        } else if (!is.null(samples) & is.null(ID)) {
            curr <- paste(samples[i], sep="")
        }
        name_vector[i] <- curr
    }
    names(final) <- name_vector
    for (i in seq_along(final)){
      final[[i]]<-final[[i]][!duplicated(final[[i]]$barcode),]
      final[[i]]<-final[[i]][rowSums(is.na(final[[i]])) < 10, ]
      final[[i]][final[[i]] == "NA"] <- NA
    }
    if (removeNA) {
      final <- .removingNA(final)
    }
    if (removeMulti) {
      final <- .removingMulti(final)
    }
    #Adding list element names to output if samples NULL
    if(is.null(samples)) {
      names(final) <- paste0("S", seq_len(length(final)))
    }
    final
}

#' Combining the list of B cell receptor contigs into clones
#'
#' This function consolidates a list of BCR sequencing results to the level 
#' of the individual cell barcodes. Using the samples and ID parameters, 
#' the function will add the strings as prefixes to prevent issues with 
#' repeated barcodes. The resulting new barcodes will need to match the 
#' Seurat or SCE object in order to use, \code{\link{combineExpression}}. 
#' Unlike \code{\link{combineTCR}}, combineBCR produces a column 
#' \strong{CTstrict} of an index of nucleotide sequence and the 
#' corresponding V gene. This index automatically calculates the 
#' Levenshtein distance between sequences with the same V gene and will
#' index sequences using a normalized Levenshtein distance with the same 
#' ID. After which, clone clusters are called using the 
#' \code{\link[igraph]{components}} function. Clones that are clustered 
#' across multiple sequences will then be labeled with "Cluster" in the 
#' CTstrict header.
#'
#' @examples
#' #Data derived from the 10x Genomics intratumoral NSCLC B cells
#' BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
#' combined <- combineBCR(BCR, 
#'                        samples = "Patient1", 
#'                        threshold = 0.85)
#' 
#' @param input.data List of filtered contig annotations or outputs from 
#' \code{\link{loadContigs}}.
#' @param samples The labels of samples
#' @param ID The additional sample labeling (optional).
#' @param call.related.clones Use the nucleotide sequence and V gene 
#' to call related clones. Default is set to TRUE. FALSE will return 
#' a CTstrict or strict clone as V gene + amino acid sequence.
#' @param threshold The normalized edit distance to consider. The higher 
#' the number the more similarity of sequence will be used for clustering.
#' @param removeNA This will remove any chain without values.
#' @param removeMulti This will remove barcodes with greater than 2 chains.
#' @param filterMulti This option will allow for the selection of the 
#' highest-expressing light and heavy chains, if not calling related clones.
#' @param filterNonproductive This option will allow for the removal of 
#' nonproductive chains if the variable exists in the contig data. Default
#' is set to TRUE to remove nonproductive contigs.
#' @import dplyr
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return List of clones for individual cell barcodes
combineBCR <- function(input.data, 
                       samples = NULL, 
                       ID = NULL, 
                       call.related.clones = TRUE,
                       threshold = 0.85,
                       removeNA = FALSE, 
                       removeMulti = FALSE,
                       filterMulti = TRUE,
                       filterNonproductive = TRUE) {
    input.data <- .checkList(input.data)
    input.data <- .checkContigs(input.data)
    out <- NULL
    final <- list()
    chain1 <- "heavy"
    chain2 <- "light"
    for (i in seq_along(input.data)) {
        input.data[[i]] <- subset(input.data[[i]], chain %in% c("IGH", "IGK", "IGL"))
        input.data[[i]]$ID <- ID[i]
        if(c("productive") %in% colnames(input.data[[i]]) & filterNonproductive) {
          input.data[[i]] <- subset(input.data[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        }
        if (filterMulti) {
                    # Keep IGH / IGK / IGL info in save_chain
                    input.data[[i]]$save_chain <- input.data[[i]]$chain
                    # Collapse IGK and IGL chains
                    input.data[[i]]$chain <- ifelse(input.data[[i]]$chain=="IGH","IGH","IGLC")
                    input.data[[i]] <- .filteringMulti(input.data[[i]])
                    # Get back IGK / IGL distinction
                    input.data[[i]]$chain <- input.data[[i]]$save_chain
                    input.data[[i]]$save_chain <- NULL
        }
    }
    if (!is.null(samples)) {
        out <- .modifyBarcodes(input.data, samples, ID)
    } else {
        out <- input.data
    }
    for (i in seq_along(out)) { 
        data2 <- data.frame(out[[i]])
        data2 <- .makeGenes(cellType = "B", data2)
        unique_df <- unique(data2$barcode)
        Con.df <- data.frame(matrix(NA, length(unique_df), 9))
        colnames(Con.df) <- c("barcode", heavy_lines, light_lines)
        Con.df$barcode <- unique_df
        Con.df <- .parseBCR(Con.df, unique_df, data2)
        Con.df <- .assignCT(cellType = "B", Con.df)
        data3<-Con.df %>% mutate(length1 = nchar(cdr3_nt1)) %>%
            mutate(length2 = nchar(cdr3_nt2))
        final[[i]] <- data3 
    }
    dictionary <- bind_rows(final)
    if(call.related.clones) {
      IGH <- .lvCompare(dictionary, "IGH", "cdr3_nt1", threshold)
      IGLC <- .lvCompare(dictionary, "IGLC", "cdr3_nt2", threshold)
    } 
    for(i in seq_along(final)) {
      if(call.related.clones) {
        final[[i]]<-merge(final[[i]],IGH,by.x="cdr3_nt1",by.y="clone",all.x=TRUE)
        final[[i]]<-merge(final[[i]],IGLC,by.x="cdr3_nt2",by.y="clone",all.x=TRUE)
        num <- ncol(final[[i]])
        final[[i]][,"CTstrict"] <- paste0(final[[i]][,num-1],".",
              final[[i]][,"vgene1"],"_",final[[i]][,num],".",final[[i]][,"vgene2"])
      } else {
        final[[i]][,"CTstrict"] <- paste0(final[[i]][,"vgene1"], ".", final[[i]][,"cdr3_aa1"], "_", final[[i]][,"vgene2"], ".", final[[i]][,"cdr3_aa2"])
      }
        final[[i]]$sample <- samples[i]
        final[[i]]$ID <- ID[i]
        final[[i]][final[[i]] == "NA_NA" | final[[i]] == "NA;NA_NA;NA"] <- NA 
        if (!is.null(sample) & !is.null(ID)) {
          final[[i]]<- final[[i]][, c("barcode", "sample", "ID", 
              heavy_lines[c(1,2,3)], light_lines[c(1,2,3)], CT_lines)]
        }
        else if (!is.null(sample) & is.null(ID)) {
          final[[i]]<- final[[i]][, c("barcode", "sample", 
                    heavy_lines[c(1,2,3)], light_lines[c(1,2,3)], CT_lines)]
        }
    }
    names <- NULL
    for (i in seq_along(samples)) { 
        if (!is.null(samples) & !is.null(ID)) {
            c <- paste(samples[i], "_", ID[i], sep="")
        } else if (!is.null(samples) & is.null(ID)) {
            c <- paste(samples[i], sep="")
        }
        names <- c(names, c)
    }
    names(final) <- names
    for (i in seq_along(final)) {
        final[[i]] <- final[[i]][!duplicated(final[[i]]$barcode),]
        final[[i]]<-final[[i]][rowSums(is.na(final[[i]])) < 10, ]
        final[[i]][final[[i]] == "NA"] <- NA
    }
    if (removeNA) { 
      final <- .removingNA(final) 
    }
    if (removeMulti) { 
      final <- .removingMulti(final) 
    }
    #Adding list element names to output if samples NULL
    if(is.null(samples)) {
      names(final) <- paste0("S", seq_len(length(final)))
    }
    return(final) 
}

