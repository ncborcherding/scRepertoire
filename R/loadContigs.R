#' Loading the contigs derived from single-cell sequencing
#'
#' This function generates a contig list and formats the data to allow for
#' function with  \code{\link{combineTCR}} or \code{\link{combineBCR}}. If
#' using data derived from filtered outputs of 10X Genomics, there is no
#' need to use this function as the data is already compatible.
#'
#' The files that this function parses includes:
#' \itemize{
#'   \item 10X =  "filtered_contig_annotations.csv"
#'   \item AIRR = "airr_rearrangement.tsv"
#'   \item BD = "Contigs_AIRR.tsv"
#'   \item Dandelion = "all_contig_dandelion.tsv"
#'   \item Immcantation = "data.tsv"
#'   \item JSON = ".json"
#'   \item ParseBio = "barcode_report.tsv"
#'   \item MiXCR = "clones.tsv"
#'   \item Omniscope = ".csv"
#'   \item TRUST4 = "barcode_report.tsv"
#'   \item WAT3R = "barcode_results.csv"
#' }
#'
#' @examples
#' TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
#' contig.list <- loadContigs(TRUST4, format = "TRUST4")
#'
#' BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
#' contig.list <- loadContigs(BD, format = "BD")
#'
#' WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
#' contig.list <- loadContigs(WAT3R, format = "WAT3R")
#'
#' @param input The directory in which contigs are located or a list with contig
#' elements
#' @param format The format of the single-cell contig, currently supporting:
#' "10X", "AIRR", "BD", "Dandelion", "JSON", "MiXCR", "ParseBio", "Omniscope",
#' "TRUST4", "WAT3R", and "Immcantation"
#' @importFrom utils read.csv read.delim
#' @importFrom rjson fromJSON
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return List of contigs for compatibility  with \code{\link{combineTCR}} or
#' \code{\link{combineBCR}}
loadContigs <- function(input, format = "10X") {

    assert_that(is.string(input) || is.list(input) || is.data.frame(input))
    assert_that(is.string(format))
    assert_that(format %in% c(
        "10X", "AIRR", "BD", "Dandelion", "JSON", "MiXCR", "ParseBio",
        "Omniscope", "TRUST4", "WAT3R", "Immcantation"
    ))

    #Loading from directory, recursively
    rawDataDfList <- if (inherits(x = input, what = "character")) {

        format.list <- list(
            "WAT3R" = "barcode_results.csv",
            "10X" =  "filtered_contig_annotations.csv",
            "AIRR" = "airr_rearrangement.tsv",
            "Dandelion" = "all_contig_dandelion.tsv",
            "Immcantation" = "_data.tsv",
            "MiXCR" = "clones.tsv",
            "JSON" = ".json",
            "TRUST4" = "barcode_report.tsv",
            "BD" = "Contigs_AIRR.tsv",
            "Omniscope" = c("_OSB.csv", "_OST.csv"),
            "ParseBio" = "barcode_report.tsv"
        )
        file.pattern <- format.list[[format]]
        contig.files <- list.files(
            input,
            paste0("*", file.pattern, "$", collapse = "|"),
            recursive = TRUE,
            full.names = TRUE
        )

        if (length(contig.files) == 0) {
            warning("No files found in the directory")
            return(list())
        }

        reader <- if (format == "json") {
            function(x) as.data.frame(fromJSON(x))
        } else if (format %in% c("10X", "WAT3R", "Omniscope")) {
            read.csv
        } else {
            read.delim
        }

        lapply(contig.files, reader)

    } else { # handle an already loaded list of dfs / 1 df
        .checkList(input)
    }

    loadFunc <- switch(format,
        "10X" = .parse10x,
        "AIRR" = .parseAIRR,
        "Dandelion" = .parseDandelion,
        "JSON" = .parseJSON,
        "MiXCR" = .parseMiXCR,
        "TRUST4" = .parseTRUST4,
        "BD" = .parseBD,
        "WAT3R"  = .parseWAT3R,
        "Omniscope" = .parseOmniscope,
        "Immcantation" = .parseImmcantation,
        "ParseBio" = .parseParse
    )

    loadFunc(rawDataDfList)
}

#' Formats TRUST4 data
#' @importFrom stringr str_split
.parseTRUST4 <- function(df) {

    for (i in seq_along(df)) {

        colnames(df[[i]])[1] <- "barcode"
        df[[i]][df[[i]] == "*"] <- NA

        if (length(which(is.na(df[[i]]$chain1))) == length(df[[i]]$chain1)) {
            chain2 <- matrix(ncol = 7, nrow = length(df[[i]]$chain1))
        } else {
            chain2 <- str_split(df[[i]]$chain1, ",", simplify = TRUE)[, seq_len(7)]
            chain2[chain2 == "*"] <- "None"
        }
        colnames(chain2) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads")
        chain2 <- data.frame(barcode = df[[i]][, 1], chain2)

        if (length(which(is.na(df[[i]]$chain2))) == length(df[[i]]$chain2)) {
            chain1 <- matrix(ncol = 7, nrow = length(df[[i]]$chain2))
        } else {
            chain1 <- str_split(df[[i]]$chain2, ",", simplify = TRUE)[, seq_len(7)]
            chain1[chain1 == "*"] <- "None"
        }
        colnames(chain1) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads") # issue 429
        chain1 <- data.frame(barcode = df[[i]][, 1], chain1)
        data2 <- rbind(chain1, chain2)
        data2[data2 == ""] <- NA
        df[[i]] <- data2
    }

    .chain.parser(df)
}

#Grabs the chain info from v_gene
.chain.parser <- function(df) {
    lapply(df, function(x) {
        x$chain <- substr(x$v_gene, 1, 3)
        x
    })
}

#Formats wat3r data
#' @author Kyle Romine, Nick Borcherding
.parseWAT3R <- function(df) {
    for (i in seq_along(df)) {
        df[[i]][df[[i]] == ""] <- NA
        chain2 <- df[[i]][,c("BC","TRBV","TRBD","TRBJ","TRB_CDR3nuc","TRB_CDR3","TRB_nReads","TRB_CDR3_UMIcount")]
        chain2 <- data.frame(chain2[,1], chain = "TRB", chain2[,2:4], c_gene = NA, chain2[,5:8])
        colnames(chain2) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
       
        #TRA Chain 1
        chain1 <-  df[[i]][,c("BC","TRAV","TRAJ","TRA_CDR3nuc","TRA_CDR3","TRA_nReads","TRA_CDR3_UMIcount")]
        chain1 <- data.frame(chain1[,1], chain = "TRA",chain1[,2], d_gene = NA, chain1[,3], c_gene = NA, chain1[,4:7])
        colnames(chain1) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
        data2 <- rbind(chain1, chain2)
        data2[data2 == ""] <- NA
       
        #TRA Chain 2
        chain3 <-  df[[i]][,c("BC","TRAV.2","TRAJ.2","TRA.2_CDR3nuc","TRA.2_CDR3","TRA.2_nReads","TRA.2_CDR3_UMIcount")]
        chain3 <- data.frame(chain3[,1], chain = "TRA",chain3[,2],  d_gene = NA, chain3[,3], c_gene = NA, chain3[,4:7])
        colnames(chain3) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")
        data2 <- rbind(chain1, chain2, chain3)
        data2[data2 == ""] <- NA
        df[[i]] <- data2
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
       
    }
    return(df)
}

#Formats AIRR data
.parseAIRR <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- df[[i]][,c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
        colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}

#Loads 10x data
.parse10x <- function(df) {
    for (i in seq_along(df)) {
        df[[i]] <- subset(df[[i]], chain != "Multi")
        df[[i]] <- subset(df[[i]], productive %in% c(TRUE, "TRUE", "True", "true"))
        if (nrow(df[[i]]) == 0) { stop(
            "There are 0 contigs after internal filtering -
            check the contig list to see if any issues exist
            for productive chains", call. = FALSE) }
        df[[i]] <- subset(df[[i]], cdr3 != "None")
        df[[i]][df[[i]] == ""] <- NA
        df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
    return(df)
}
#Loads BD AIRR
.parseBD <- function(df) {
  for (i in seq_along(df)) {
    df[[i]] <- df[[i]][,c("cell_id","locus","v_call","d_call","j_call", "c_call", "cdr3","cdr3_aa","consensus_count", "productive")]
    colnames(df[[i]]) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "productive")
    df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
  }
  return(df)
}

.parseOmniscope <- function(df) {
  for (i in seq_along(df)) {
    if("c_call" %in% colnames(df[[i]])) {
      df[[i]] <- df[[i]][,c("cell_id", "locus", "umi_count", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "productive")
    } else { #TCR contigs do not include C gene
      df[[i]] <- df[[i]][,c("cell_id", "locus", "umi_count", "v_call", "d_call", "j_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "cdr3_nt", "cdr3", "productive")
      df[[i]][,"c_gene"] <- NA
    df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
    }
  }
  return(df)
}

.parseJSON <- function(df) {
  for (i in seq_along(df)) {
    df[[i]] <- do.call(rbind, df[[i]])
    df[[i]][df[[i]] == ""] <- NA
    df[[i]] <- as.data.frame(df[[i]])
    df[[i]] <- df[[i]][,c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
    colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
  }
  return(df)
}

.parseMiXCR <- function(df) {
  for (i in seq_along(df)) {
    df[[i]][df[[i]] == ""] <- NA
    df[[i]] <- as.data.frame(df[[i]])
    df[[i]] <- df[[i]][,c("tagValueCELL", "topChains", "readCount", "allVHitsWithScore",   "allDHitsWithScore",   "allJHitsWithScore",  "allCHitsWithScore", "nSeqCDR3", "aaSeqCDR3")]
    colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
  }
  return(df)
}

#' @importFrom stringr str_split
.parseImmcantation<- function(df) {
  for (i in seq_along(df)) {
    df[[i]][df[[i]] == ""] <- NA
    df[[i]] <- as.data.frame(df[[i]])
    if("c_call" %in% colnames(df[[i]])) {
      df[[i]] <- df[[i]][,c("sequence_id", "locus", "consensus_count",  "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "productive")
    } else {
      df[[i]] <- df[[i]][,c("sequence_id", "locus", "consensus_count",  "v_call", "d_call", "j_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "cdr3_nt", "cdr3", "productive")
      df[[i]][,"c_gene"] <- NA
    }
    df[[i]]$barcode <- str_split(df[[i]][,"barcode"], "_", simplify = TRUE)[,1]
  }
  return(df)
}

.parseParse <- function(df) {
  for (i in seq_along(df)) {
    df[[i]][df[[i]] == ""] <- NA
    df[[i]][df[[i]] == "NaN"] <- NA
    df[[i]][df[[i]] == "nan"] <- NA
    df[[i]] <- as.data.frame(df[[i]])
    TRA.1 <- df[[i]][,c("Barcode", "TRA_V", "TRA_D", "TRA_J", "TRA_C", "TRA_cdr3_aa", "TRA_read_count", "TRA_transcript_count")]
    TRA.2 <- df[[i]][,c("Barcode", "secondary_TRA_V", "secondary_TRA_D", "secondary_TRA_J", "secondary_TRA_C", "secondary_TRA_cdr3_aa", "secondary_TRA_read_count", "secondary_TRA_transcript_count")]
    colnames(TRA.1) <- 1:8
    colnames(TRA.2) <- 1:8
    TRA <- rbind(TRA.1, TRA.2)
    TRA$chain <- "TRA"
   
    TRB.1 <- df[[i]][,c("Barcode", "TRB_V", "TRB_D", "TRB_J", "TRB_C", "TRB_cdr3_aa", "TRB_read_count", "TRB_transcript_count")]
    TRB.2 <- df[[i]][,c("Barcode", "secondary_TRB_V", "secondary_TRB_D", "secondary_TRB_J", "secondary_TRB_C", "secondary_TRB_cdr3_aa", "secondary_TRB_read_count", "secondary_TRB_transcript_count")]
    colnames(TRB.1) <- 1:8
    colnames(TRB.2) <- 1:8
    TRB <- rbind(TRB.1, TRB.2)
    TRB$chain <- "TRB"
   
    data2 <- rbind(TRA, TRB)
    data2 <- data2[rowSums(is.na(data2[2:8])) != 7, ]
    colnames(data2) <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "reads", "umis", "chain")
    data2$cdr3_nt <- NA
    data2 <- data2[,c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")]
   
    df[[i]] <- data2
    df[[i]] <- df[[i]][with(df[[i]], order(reads, chain)),]
  }
  return(df)
}

.parseDandelion <- function(df) {
  for (i in seq_along(df)) {
    df[[i]] <- df[[i]][,c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "productive")]
    colnames(df[[i]]) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "productive")
  }
  return(df)
}