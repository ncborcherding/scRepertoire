#' Load Immune Receptor Sequencing Contigs
#'
#' @description
#' This function loads and processes contig data from various single-cell 
#' immune receptor sequencing formats. It reads data from a directory 
#' (recursively) or from an already loaded list/data frame, transforms it to a 
#' common structure, and returns a list of contigs ready for downstream analysis 
#' with [combineTCR()] or [combineBCR()].
#'
#' Supported file formats and their expected file names:
#' - `10X`: "filtered_contig_annotations.csv"
#' - `AIRR`: "airr_rearrangement.tsv"
#' - `BD`: "Contigs_AIRR.tsv"
#' - `Dandelion`: "all_contig_dandelion.tsv"
#' - `Immcantation`: "_data.tsv" (or similar)
#' - ``JSON`: ".json"
#' - `ParseBio`: "barcode_report.tsv"
#' - `MiXCR`: "clones.tsv"
#' - `TRUST4`: "barcode_report.tsv"
#' - `WAT3R`: "barcode_results.csv"
#'
#' @param input A directory path containing contig files or a list/data frame 
#' of pre-loaded contig data.
#' @param format A string specifying the data format. Must be one of:
#' `auto`, `10X`, `AIRR`, `BD`, `Dandelion`, `JSON`, `MiXCR`, `ParseBio`, 
#' `TRUST4`, `WAT3R`, or `Immcantation`. If "auto", the function attempts 
#' automatic format detection.
#'
#' @return A list of contigs formatted for use with [combineTCR()] or 
#' [combineBCR()]. Rows containing only NA values (aside from the barcode) 
#' are dropped.
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
#' @importFrom utils read.csv read.delim
#' @importFrom rjson fromJSON
#' @export
loadContigs <- function(input, 
                        format = "10X") {
  valid_formats <- c("auto", "10X", "AIRR", "BD", "Dandelion", "JSON",
                     "MiXCR", "ParseBio", "TRUST4", "WAT3R", "Immcantation")
  
  # Define mapping from formats to file patterns (for directory inputs)
  format_patterns <- list(
    "WAT3R"       = "barcode_results.csv",
    "10X"         = "filtered_contig_annotations.csv",
    "AIRR"        = "airr_rearrangement.tsv",
    "Dandelion"   = "all_contig_dandelion.tsv",
    "Immcantation"= "_data.tsv",
    "MiXCR"       = "clones.tsv",
    "JSON"        = ".json",
    "TRUST4"      = "barcode_report.tsv",
    "BD"          = "Contigs_AIRR.tsv",
    "ParseBio"    = "barcode_report.tsv"
  )
  
  # Determine data input type and detect format if necessary
  if (is.character(input)) {
    # Input is a directory; use file name patterns for detection if "auto"
    if (format == "auto") {
      candidates <- c()
      for (fmt in names(format_patterns)) {
        pattern <- format_patterns[[fmt]]
        pattern_regex <- paste0(".*(", pattern, ")$")
        files_found <- list.files(
          path = input,
          pattern = pattern_regex,
          recursive = TRUE,
          full.names = TRUE
        )
        if (length(files_found) > 0) {
          candidates <- c(candidates, fmt)
        }
      }
      if (length(candidates) == 0) {
        warning("No contig files found in the specified directory using any known patterns.")
        return(list())
      } else if (length(unique(candidates)) > 1) {
        warning("Ambiguous format detection. Multiple candidate formats found: ",
                paste(unique(candidates), collapse = ", "),
                ". Defaulting to first candidate: ", unique(candidates)[1])
        detected_format <- unique(candidates)[1]
      } else {
        detected_format <- unique(candidates)[1]
        message("Automatically detected format: ", detected_format)
      }
    } else {
      detected_format <- format
    }
    
    # Use the detected or specified format to list files
    file_pattern <- format_patterns[[detected_format]]
    pattern_regex <- paste0(".*(", file_pattern, ")$")
    contig_files <- list.files(
      path = input,
      pattern = pattern_regex,
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(contig_files) == 0) {
      warning("No contig files found in the specified directory with pattern ", file_pattern, ".")
      return(list())
    }
    
    # Choose a reader based on format
    reader <- switch(detected_format,
                     "JSON" = function(x) as.data.frame(fromJSON(x)),
                     "10X" = read.csv,
                     "WAT3R" = read.csv,
                     read.delim)
    rawDataDfList <- lapply(contig_files, reader)
  } else {
    # Input is pre-loaded data; ensure it is a list of data frames
    rawDataDfList <- .checkList(input)
    if (format == "auto") {
      detected_format <- .detectFormatFromData(rawDataDfList)
      message("Automatically detected format from data: ", detected_format)
    } else {
      detected_format <- format
    }
  }
  
  # Select the appropriate parser function based on the detected format
  loadFunc <- switch(detected_format,
                     "10X" = .parse10x,
                     "AIRR" = .parseAIRR,
                     "Dandelion" = .parseDandelion,
                     "JSON" = .parseJSON,
                     "MiXCR" = .parseMiXCR,
                     "TRUST4" = .parseTRUST4,
                     "BD" = .parseBD,
                     "WAT3R" = .parseWAT3R,
                     "Immcantation" = .parseImmcantation,
                     "ParseBio" = .parseParse)
  
  # Process and clean contig data, dropping rows that are entirely NA (except for barcode)
  cleaned_list <- loadFunc(rawDataDfList)
  .rmAllNaRowsFromLoadContigs(cleaned_list)
}

# Auto-detect format from a pre-loaded data object
.detectFormatFromData <- function(dfList) {
  # Use the first data frame for detection
  df <- if (is.data.frame(dfList)) dfList else dfList[[1]]
  cn <- colnames(df)
  
  if (all(c("chain", "productive", "cdr3") %in% cn)) {
    return("10X")
  } else if (all(c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa") %in% cn)) {
    return("AIRR")
  } else if (all(c("cell_id", "locus", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "consensus_count", "productive") %in% cn)) {
    # Could be BD or Dandelion – defaulting to BD here.
    return("BD")
  } else if (all(c("tagValueCELL", "topChains", "readCount", "allVHitsWithScore") %in% cn)) {
    return("MiXCR")
  } else if (all(c("barcode", "chain1", "chain2") %in% cn)) {
    return("TRUST4")
  } else if (all(c("BC", "TRBV", "TRBD", "TRBJ") %in% cn)) {
    return("WAT3R")
  } else if (all(c("sequence_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "cdr3") %in% cn)) {
    if ("productive" %in% cn) return("Immcantation")
  } else if (all(c("Barcode", "TRA_V") %in% cn) || all(c("Barcode", "IGH_V") %in% cn)) {
    return("ParseBio")
  }
  stop("Could not auto-detect format from the provided data structure.")
}

# Remove rows with all NA values (ignoring the "barcode" column)
.rmAllNaRowsFromLoadContigs <- function(dfList) {
  cols <- setdiff(colnames(dfList[[1]]), "barcode")
  lapply(dfList, function(df) {
    df[rowSums(!is.na(df[, cols, drop = FALSE])) > 0, , drop = FALSE]
  })
}

# Helper function to replace empty strings with NA
.sanitize_empty <- function(df) {
  df[df == ""] <- NA
  df
}

# Helper function to order a data frame by "reads" and "chain" (if present)
.order_df <- function(df) {
  if (all(c("reads", "chain") %in% colnames(df))) {
    df[order(df$reads, df$chain), ]
  } else {
    df
  }
}

## --- Parsing Functions ---
# (Parsing functions below remain unchanged for completeness)

.parseTRUST4 <- function(df_list) {
  split_and_pad <- function(x, n = NULL) {
    parts <- unlist(strsplit(x, ","))
    if (is.null(n)) return(parts)
    if (length(parts) < n) {
      return(c(parts, rep("", n - length(parts))))
    } else if (length(parts) > n) {
      return(parts[seq_len(n)])
    } else {
      return(parts)
    }
  }
  
  processChain <- function(data, chain_col) {
    if (all(is.na(data[[chain_col]]))) {
      chain <- matrix(NA, ncol = 7, nrow = length(data[[chain_col]]))
    } else {
      chain <- t(vapply(
        data[[chain_col]], 
        split_and_pad, 
        FUN.VALUE = character(7), 
        n = 7,
        USE.NAMES = FALSE
      ))
      chain[chain == "*"] <- "None"
    }
    colnames(chain) <- c("v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads")
    data.frame(barcode = data$barcode, chain, stringsAsFactors = FALSE)
  }
  
  formatted <- lapply(df_list, function(data) {
    colnames(data)[1] <- "barcode"
    data <- .sanitize_empty(data)
    chain1 <- processChain(data, "chain2")
    chain2 <- processChain(data, "chain1")
    combined_data <- rbind(chain1, chain2)
    .sanitize_empty(combined_data)
  })
  
  lapply(formatted, function(x) {
    x$chain <- substr(x$v_gene, 1, 3)
    x
  })
}

.parseWAT3R <- function(df_list) {
  lapply(df_list, function(df) {
    df <- .sanitize_empty(df)
    chain2 <- data.frame(
      barcode = df$BC,
      chain = "TRB",
      v_gene = df$TRBV,
      d_gene = df$TRBD,
      j_gene = df$TRBJ,
      c_gene = NA,
      cdr3_nt = df$TRB_CDR3nuc,
      cdr3 = df$TRB_CDR3,
      reads = df$TRB_nReads,
      umis = df$TRB_CDR3_UMIcount,
      stringsAsFactors = FALSE
    )
    chain1 <- data.frame(
      barcode = df$BC,
      chain = "TRA",
      v_gene = df$TRAV,
      d_gene = NA,
      j_gene = df$TRAJ,
      c_gene = NA,
      cdr3_nt = df$TRA_CDR3nuc,
      cdr3 = df$TRA_CDR3,
      reads = df$TRA_nReads,
      umis = df$TRA_CDR3_UMIcount,
      stringsAsFactors = FALSE
    )
    chain3 <- data.frame(
      barcode = df$BC,
      chain = "TRA",
      v_gene = df[["TRAV.2"]],
      d_gene = NA,
      j_gene = df[["TRAJ.2"]],
      c_gene = NA,
      cdr3_nt = df[["TRA.2_CDR3nuc"]],
      cdr3 = df[["TRA.2_CDR3"]],
      reads = df[["TRA.2_nReads"]],
      umis = df[["TRA.2_CDR3_UMIcount"]],
      stringsAsFactors = FALSE
    )
    combined <- rbind(chain1, chain2, chain3)
    .order_df(.sanitize_empty(combined))
  })
}

.parseAIRR <- function(df_list) {
  lapply(df_list, function(df) {
    df <- df[, c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
    colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
    .order_df(df)
  })
}

.parse10x <- function(df_list) {
  lapply(df_list, function(df) {
    df <- subset(df, chain != "Multi")
    df <- subset(df, productive %in% c(TRUE, "TRUE", "True", "true"))
    if (nrow(df) == 0) {
      stop("There are 0 contigs after filtering. Check the contig list for issues with productive chains.", call. = FALSE)
    }
    df <- subset(df, cdr3 != "None")
    df <- .sanitize_empty(df)
    .order_df(df)
  })
}

.parseBD <- function(df_list) {
  lapply(df_list, function(df) {
    df <- df[, c("cell_id", "locus", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "consensus_count", "productive")]
    colnames(df) <- c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "productive")
    .order_df(df)
  })
}

.parseJSON <- function(df_list) {
  lapply(df_list, function(df) {
    df <- do.call(rbind, df)
    df <- .sanitize_empty(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- df[, c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")]
    colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
    df
  })
}

.parseMiXCR <- function(df_list) {
  lapply(df_list, function(df) {
    df <- .sanitize_empty(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- df[, c("tagValueCELL", "topChains", "readCount", "allVHitsWithScore",
                 "allDHitsWithScore", "allJHitsWithScore", "allCHitsWithScore",
                 "nSeqCDR3", "aaSeqCDR3")]
    colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene",
                      "j_gene", "c_gene", "cdr3_nt", "cdr3")
    df
  })
}

.parseImmcantation <- function(df_list) {
  lapply(df_list, function(df) {
    df <- .sanitize_empty(df)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    if ("c_call" %in% colnames(df)) {
      df <- df[, c("sequence_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "productive")
    } else {
      df <- df[, c("sequence_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "cdr3", "cdr3_aa", "productive")]
      colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "cdr3_nt", "cdr3", "productive")
      df$c_gene <- NA
    }
    df$barcode <- vapply(strsplit(df$barcode, "_"), `[`, 1, FUN.VALUE = character(1))
    df
  })
}

.parseParse <- function(df_list) {
  lapply(df_list, function(df) {
    df <- .sanitize_empty(df)
    df[df %in% c("NaN", "nan")] <- NA
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    
    isTcell <- any(c("TRA_V", "TRA_D", "TRA_J") %in% colnames(df))
    if (isTcell) {
      TRA1 <- df[, c("Barcode", "TRA_V", "TRA_D", "TRA_J", "TRA_C", "TRA_cdr3_aa", "TRA_read_count", "TRA_transcript_count")]
      TRA2 <- df[, c("Barcode", "secondary_TRA_V", "secondary_TRA_D", "secondary_TRA_J", "secondary_TRA_C", "secondary_TRA_cdr3_aa", "secondary_TRA_read_count", "secondary_TRA_transcript_count")]
      colnames(TRA1) <- colnames(TRA2) <- as.character(1:8)
      TRA <- rbind(TRA1, TRA2)
      TRA$chain <- "TRA"
      
      TRB1 <- df[, c("Barcode", "TRB_V", "TRB_D", "TRB_J", "TRB_C", "TRB_cdr3_aa", "TRB_read_count", "TRB_transcript_count")]
      TRB2 <- df[, c("Barcode", "secondary_TRB_V", "secondary_TRB_D", "secondary_TRB_J", "secondary_TRB_C", "secondary_TRB_cdr3_aa", "secondary_TRB_read_count", "secondary_TRB_transcript_count")]
      colnames(TRB1) <- colnames(TRB2) <- as.character(1:8)
      TRB <- rbind(TRB1, TRB2)
      TRB$chain <- "TRB"
      
      combined <- rbind(TRA, TRB)
    } else {
      IGH1 <- df[, c("Barcode", "IGH_V", "IGH_D", "IGH_J", "IGH_C", "IGH_cdr3_aa", "IGH_read_count", "IGH_transcript_count")]
      IGH2 <- df[, c("Barcode", "secondary_IGH_V", "secondary_IGH_D", "secondary_IGH_J", "secondary_IGH_C", "secondary_IGH_cdr3_aa", "secondary_IGH_read_count", "secondary_IGH_transcript_count")]
      colnames(IGH1) <- colnames(IGH2) <- as.character(1:8)
      IGH <- rbind(IGH1, IGH2)
      IGH$chain <- "IGH"
      
      IGK1 <- df[, c("Barcode", "IGK_V", "IGK_D", "IGK_J", "IGK_C", "IGK_cdr3_aa", "IGK_read_count", "IGK_transcript_count")]
      IGK2 <- df[, c("Barcode", "secondary_IGK_V", "secondary_IGK_D", "secondary_IGK_J", "secondary_IGK_C", "secondary_IGK_cdr3_aa", "secondary_IGK_read_count", "secondary_IGK_transcript_count")]
      colnames(IGK1) <- colnames(IGK2) <- as.character(1:8)
      IGK <- rbind(IGK1, IGK2)
      IGK$chain <- "IGK"
      
      IGL1 <- df[, c("Barcode", "IGL_V", "IGL_D", "IGL_J", "IGL_C", "IGL_cdr3_aa", "IGL_read_count", "IGL_transcript_count")]
      IGL2 <- df[, c("Barcode", "secondary_IGL_V", "secondary_IGL_D", "secondary_IGL_J", "secondary_IGL_C", "secondary_IGL_cdr3_aa", "secondary_IGL_read_count", "secondary_IGL_transcript_count")]
      colnames(IGL1) <- colnames(IGL2) <- as.character(1:8)
      IGL <- rbind(IGL1, IGL2)
      IGL$chain <- "IGL"
      
      combined <- rbind(IGH, IGK, IGL)
    }
    combined <- combined[rowSums(is.na(combined[, 2:8])) != 7, ]
    colnames(combined) <- c("barcode", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3", "reads", "umis", "chain")
    combined$cdr3_nt <- NA
    combined <- combined[, c("barcode", "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "reads", "umis")]
    .order_df(combined)
  })
}

.parseDandelion <- function(df_list) {
  lapply(df_list, function(df) {
    df <- df[, c("cell_id", "locus", "consensus_count", "v_call", "d_call", "j_call", "c_call", "cdr3", "cdr3_aa", "productive")]
    colnames(df) <- c("barcode", "chain", "reads", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3", "productive")
    df
  })
}
