#' Export Clonal Data in Various Formats
#'
#' Exports clonal information (gene sequences, amino acids, nucleotides) from
#' scRepertoire objects into a file or a data frame. The output format can be
#' tailored for compatibility with different analysis workflows.
#'
#' @details
#' The `format` parameter determines the structure of the output:
#' \itemize{
#'   \item{`paired`}: Exports a data frame where each row represents a barcode,
#'     with paired chain information (amino acid, nucleotide, genes) in separate
#'     columns.
#'   \item{`airr`}: Exports a data frame that adheres to the Adaptive Immune
#'     Receptor Repertoire (AIRR) Community format, with each row representing
#'     a single receptor chain.
#'   \item{`TCRMatch`}: Exports a data frame specifically for the TCRMatch
#'     algorithm, containing the TRB chain amino acid sequence and clonal
#'     frequency.
#'   \item{`tcrpheno`}: Exports a data frame compatible with the `tcrpheno`
#'     pipeline, with TRA and TRB chains in separate columns.
#'   \item{`immunarch`}: Exports a list containing a data frame and metadata
#'     formatted for use with the `immunarch` package.
#' }
#'
#' @param input.data The product of `combineTCR()`, `combineBCR()`, or
#'   `combineExpression()`.
#' @param format The format for exporting clones.
#'   Options are: `paired`, `airr`, `TCRMatch`, `tcrpheno`, `immunarch`.
#' @param group.by The variable in the metadata to use for grouping. If `NULL`,
#'   data will be grouped by the sample names.
#' @param write.file If `TRUE` (default), saves the output to a CSV file. If
#'   `FALSE`, returns the data frame or list to the R environment.
#' @param dir The directory where the output file will be saved.
#'   Defaults to the current working directory.
#' @param file.name The name of the file to be saved.
#'
#' @return A data frame or list in the specified format, either returned to the
#'   R environment or saved as a CSV file.
#'
#' @author Jonathan Noonan, Nick Borcherding
#'
#' @importFrom utils write.csv
#' @export
#' @concept Loading_and_Processing_Contigs
#'
#' @examples
#' \dontrun{
#' #Making combined contig data
#' combined <- combineTCR(contig_list,
#'                        samples = c("P17B", "P17L", "P18B", "P18L",
#'                                    "P19B", "P19L", "P20B", "P20L"))
#'
#' # Export as a paired data frame and save to a file
#' exportClones(combined, format = "paired", file.name = "paired_clones.csv")
#'
#' # Return an AIRR-formatted data frame to the environment
#' airr_df <- exportClones(combined, format = "airr", write.file = FALSE)
#' }
exportClones <- function(input.data,
                         format = "paired",
                         group.by = NULL,
                         write.file = TRUE,
                         dir = NULL,
                         file.name = "clones.csv") {
  
  # Validate format parameter
  format <- match.arg(format, c("paired", "airr", "TCRMatch", "tcrpheno", "immunarch"))
  
  # Select the appropriate internal export function
  exportFunc <- switch(format,
                       "paired"     = .pairedExport,
                       "airr"       = .airrExport,
                       "TCRMatch"   = .tcrMatchExport,
                       "tcrpheno"   = .tcrPhenoExport,
                       "immunarch"  = .immunarchExport
  )
  
  # Generate the data matrix/list
  output_data <- exportFunc(input.data, group.by)
  
  # Replace string "NA" with actual NA values
  if (is.data.frame(output_data)) {
    output_data[output_data == "NA"] <- NA
  }
  
  if (!write.file) {
    return(output_data)
  }
  
  # Handle file writing
  if (is.null(dir)) {
    dir <- "."
  }
  filepath <- file.path(dir, file.name)
  
  # Immunarch format returns a list of data frames. To save as a single CSV,
  # we bind them together and add a 'Sample' identifier column.
  if (format == "immunarch") {
    bound_data <- dplyr::bind_rows(output_data$data, .id = "Sample")
    write.csv(bound_data, file = filepath, row.names = FALSE)
  } else {
    write.csv(output_data, file = filepath)
  }
}


#' Split String and Pad to a Fixed-Width Matrix
#' @keywords internal
._split_and_pad <- function(x, split, n_cols) {
  s <- strsplit(x, split)
  # Create a matrix by safely subsetting each list element up to n_cols
  mat <- t(sapply(s, `[`, seq_len(n_cols)))
  mat[mat == "NA"] <- NA
  return(mat)
}

#' Bind a List of Contig Data Frames and Add Grouping Variable
#' @keywords internal
._bind_contig_list <- function(contig_list) {
  list_names <- names(contig_list)
  if (is.null(list_names)) {
    list_names <- as.character(seq_along(contig_list))
  }
  # Efficiently add group names before binding
  contig_list <- mapply(function(df, name) transform(df, group = name),
                        contig_list, list_names, SIMPLIFY = FALSE)
  # Combine all data frames
  return(do.call(rbind, contig_list))
}

# remap a genes matrix (n x 4) to V/D/J/C with NA for missing D in 3-gene chains
.remap_VDJC <- function(genes_mat) {
  # Coerce to a 4-col character matrix, pad if needed
  if (is.null(genes_mat)) stop(".remap_VDJC: genes_mat is NULL")
  genes_mat <- as.matrix(genes_mat)
  storage.mode(genes_mat) <- "character"
  if (ncol(genes_mat) < 4) {
    pad <- matrix(NA_character_, nrow(genes_mat), 4 - ncol(genes_mat))
    genes_mat <- cbind(genes_mat, pad)
  }
  
  # Clean tokens: drop everything after first ';', trim, map "", "NA" -> NA
  clean_token <- function(x) {
    x <- sub(";.*", "", x, perl = TRUE)
    x <- trimws(x)
    x[x == "" | x == "NA"] <- NA_character_
    x
  }
  genes_mat[] <- clean_token(genes_mat)
  
  vcol <- genes_mat[, 1]
  
  # Detect 3-gene layouts (TRA/TRG/IGK/IGL); handle NAs
  layout3 <- !is.na(vcol) & grepl("^(TRAV|TRGV|IGKV|IGLV)", vcol, perl = TRUE)
  
  # Initialize as 4-gene (V,D,J,C)
  V <- vcol
  D <- genes_mat[, 2]
  J <- genes_mat[, 3]
  C <- genes_mat[, 4]
  
  # Overwrite 3-gene rows: V, (no D), J <- col2, C <- col3
  if (any(layout3, na.rm = TRUE)) {
    idx <- which(layout3)
    D[idx] <- NA_character_
    J[idx] <- genes_mat[idx, 2]
    C[idx] <- genes_mat[idx, 3]
  }
  
  # Ensure character vectors with row length
  stopifnot(length(V) == nrow(genes_mat),
            length(D) == nrow(genes_mat),
            length(J) == nrow(genes_mat),
            length(C) == nrow(genes_mat))
  
  list(V = V, D = D, J = J, C = C)
}



# Format-Specific Export Functions 
.tcrMatchExport <- function(input.data, group.by) {
  # Wrangle data to get TRB chain information
  input.data.list <- .dataWrangle(input.data, group.by, "CTgene", "TRB")
  
  # Process each list element to ensure proper chain format
  input.data.list <- lapply(input.data.list, function(df) {
    df <- .offTheChain(df, "TRB", "CTaa", check = FALSE)
    .offTheChain(df, "TRB", "CTnt", check = FALSE)
  })
  
  bound_data <- ._bind_contig_list(input.data.list)
  # Ensure only TRB chains are included
  bound_data <- bound_data[substr(bound_data$CTgene, 1, 3) == "TRB", ]
  
  # Format for TCRMatch: remove leading character from CDR3 sequence
  mat <- data.frame(
    chain2_aa = substring(bound_data[, "CTaa"], 2),
    group = bound_data[, "group"]
  )
  mat <- na.omit(mat)
  
  # Calculate clonal frequency within each group
  mat$clonalFrequency <- ave(rep(1, nrow(mat)), mat[, c("group", "chain2_aa")], FUN = sum)
  mat <- unique(mat)
  
  return(mat)
}

#' @noRd
.pairedExport <- function(input.data, group.by) {
  input.data.list <- .dataWrangle(input.data, group.by, "CTgene", "both")
  bound_data <- ._bind_contig_list(input.data.list)
  
  # Split concatenated strings into separate columns
  genes <- ._split_and_pad(bound_data[, "CTgene"], "_", 2)
  aa    <- ._split_and_pad(bound_data[, "CTaa"], "_", 2)
  nt    <- ._split_and_pad(bound_data[, "CTnt"], "_", 2)
  
  mat <- data.frame(
    row.names    = bound_data[, "barcode"],
    chain1_aa    = aa[, 1],
    chain1_nt    = nt[, 1],
    chain1_genes = genes[, 1],
    chain2_aa    = aa[, 2],
    chain2_nt    = nt[, 2],
    chain2_genes = genes[, 2],
    group        = bound_data[, "group"]
  )
  return(mat)
}

#' @noRd
.tcrPhenoExport <- function(input.data, group.by) {
  input.data.list <- .dataWrangle(input.data, NULL, "CTgene", "both")
  dat <- do.call(rbind, unname(input.data.list))
  
  # Split sequence and gene columns
  aa_split   <- ._split_and_pad(dat$CTaa, "_", 2)
  nt_split   <- ._split_and_pad(dat$CTnt, "_", 2)
  gene_split <- ._split_and_pad(dat$CTgene, "_", 2)
  
  # Further split gene segments
  tcr1_genes <- ._split_and_pad(gene_split[, 1], "[.]", 2)
  tcr2_genes <- ._split_and_pad(gene_split[, 2], "[.]", 3)[, c(1, 3)]
  
  contigs <- data.frame(
    cell        = rownames(dat),
    TCRA_cdr3aa = aa_split[, 1],
    TCRA_vgene  = tcr1_genes[, 1],
    TCRA_jgene  = tcr1_genes[, 2],
    TCRA_cdr3nt = nt_split[, 1],
    TCRB_cdr3aa = aa_split[, 2],
    TCRB_vgene  = tcr2_genes[, 1],
    TCRB_jgene  = tcr2_genes[, 2],
    TCRB_cdr3nt = nt_split[, 2]
  )
  
  # Remove multi-chain artifacts (keep first sequence)
  columns_to_clean <- c("TCRA_cdr3aa", "TCRA_cdr3nt", "TCRB_cdr3aa", "TCRB_cdr3nt")
  contigs[columns_to_clean] <- lapply(contigs[columns_to_clean], function(x) {
    gsub(";.*", "", x)
  })
  
  # Clean up empty or "NA" strings
  contigs[contigs == "" | contigs == "NA"] <- NA
  return(contigs)
}


#' @noRd
.airrExport <- function(input.data, group.by) {
  input.data.list <- .dataWrangle(input.data, group.by, "CTgene", "both")
  bound_data <- ._bind_contig_list(input.data.list)
  
  # Helper to clean multi-chain entries by taking the first one
  clean_multichain <- function(x) {
    gsub(";.*", "", x)
  }
  
  # Split combined fields for both chains at once
  genes <- ._split_and_pad(bound_data$CTgene, "_", 2)
  aa <- ._split_and_pad(clean_multichain(bound_data$CTaa), "_", 2)
  nt <- ._split_and_pad(clean_multichain(bound_data$CTnt), "_", 2)
  
  # Process chain 1 and chain 2 gene calls
  genes1 <- ._split_and_pad(genes[, 1], "[.]", 4)
  genes2 <- ._split_and_pad(genes[, 2], "[.]", 4)
  
  # Function to create a data frame for a single chain
  format_chain_df <- function(cell_ids, locus_prefix, genes_mat, nt_seq, aa_seq) {
    # Determine if chain is heavy/beta (has D gene)
    is_heavy_or_beta <- locus_prefix %in% c("IGH", "TRB", "TRD")
    is_heavy_or_beta <- names(table(is_heavy_or_beta))[1]
    
    df <- data.frame(
      locus = paste0(locus_prefix, "1"),
      v_call = genes_mat[, 1],
      d_call = if (is_heavy_or_beta) genes_mat[, 2] else NA,
      j_call = if (is_heavy_or_beta) genes_mat[, 3] else genes_mat[, 2],
      c_call = if (is_heavy_or_beta) genes_mat[, 4] else genes_mat[, 3],
      junction = nt_seq,
      junction_aa = aa_seq
    )
    return(df)
  }
  
  # Locus prefixes for chain 1 and 2
  locus1_prefix <- substr(genes[, 1], 1, 3)
  locus2_prefix <- substr(genes[, 2], 1, 3)
  
  # Create data frames for each chain
  df1 <- format_chain_df(bound_data$barcode, locus1_prefix, genes1, nt[, 1], aa[, 1])
  df2 <- format_chain_df(bound_data$barcode, locus2_prefix, genes2, nt[, 2], aa[, 2])
  
  # Add cell_id to both
  df1$cell_id <- bound_data$barcode
  df2$cell_id <- bound_data$barcode
  
  # Combine chain data frames
  mat <- rbind(df1, df2)
  
  # Reorder columns to AIRR standard
  airr_cols <- c("cell_id", "locus", "v_call", "d_call", "j_call", "c_call", "junction", "junction_aa")
  mat <- mat[, airr_cols]
  
  # Remove rows that are entirely NA except for cell_id
  mat <- mat[rowSums(is.na(mat[, -1])) < (ncol(mat) - 1), ]
  
  return(mat)
}


#' @noRd
.immunarchExport <- function(input.data, group.by) {
  df_list <- .dataWrangle(input.data, group.by, "CTgene", "both")
  meta <- data.frame(Sample = names(df_list))
  
  data_out <- lapply(df_list, function(x) {
    # Summarize by clonotype (CTstrict)
    result <- x %>%
      dplyr::group_by(.data[["CTstrict"]]) %>%
      dplyr::summarise(
        Clones  = dplyr::n(),
        barcode = paste(barcode, collapse = ";"),
        CTaa    = dplyr::first(.data[["CTaa"]]),
        CTnt    = dplyr::first(.data[["CTnt"]]),
        CTgene  = dplyr::first(.data[["CTgene"]]),
        .groups = "drop"
      ) %>%
      dplyr::mutate(Proportion = Clones / sum(Clones))
    
    # Determine chain type from first half
    first_genes <- ._split_and_pad(result$CTgene, "_", 2)[, 1]
    chain_type  <- substr(first_genes, 1, 3)
    
    # Impute Missing Chain info (fix condition)
    if (any(is.na(chain_type)) && length(unique(chain_type)) == 2) {
      tbl.store <- rev(sort(table(chain_type, useNA = "ifany")))
      if (length(tbl.store) == 2) {
        chain_type[is.na(chain_type)] <- names(tbl.store)[1]
      }
    }
    
    # Split sequence/gene info
    aa_split   <- ._split_and_pad(result$CTaa,   "_", 2)
    nt_split   <- ._split_and_pad(result$CTnt,   "_", 2)
    gene_split <- ._split_and_pad(result$CTgene, "_", 2)
    
    .is_layout3_vec <- function(v_gene_first_token) {
      grepl("^(TRAV|TRGV|IGKV|IGLV)", v_gene_first_token %||% "", perl = TRUE)
    }
    
    # Adaptive per-row ordering: alpha/light first
    first.gene.pos <- ifelse(chain_type %in% c("TRA","TRG","IGK","IGL"), 1L, 2L)
    na_pos <- is.na(first.gene.pos)
    if (any(na_pos)) {
      left_is3  <- .is_layout3_vec(gsub("\\..*$", "", gene_split[, 1]))
      right_is3 <- .is_layout3_vec(gsub("\\..*$", "", gene_split[, 2]))
      first.gene.pos[na_pos] <- ifelse(left_is3[na_pos], 1L, 2L)
    }
    second.gene.pos <- 3L - first.gene.pos
    
    n    <- nrow(result)
    ridx <- seq_len(n)
    
    # Per-row pick halves
    aa_first   <- as.character(aa_split[cbind(ridx, first.gene.pos)])
    aa_second  <- as.character(aa_split[cbind(ridx, second.gene.pos)])
    nt_first   <- as.character(nt_split[cbind(ridx, first.gene.pos)])
    nt_second  <- as.character(nt_split[cbind(ridx, second.gene.pos)])
    g_first_s  <- as.character(gene_split[cbind(ridx, first.gene.pos)])
    g_second_s <- as.character(gene_split[cbind(ridx, second.gene.pos)])
    
    # Chain presence flags (entire half missing?)
    present1 <- !(is.na(g_first_s)  | g_first_s  == "" | g_first_s  == "NA")
    present2 <- !(is.na(g_second_s) | g_second_s == "" | g_second_s == "NA")
    
    # Tokenize per selected half
    genes1 <- ._split_and_pad(g_first_s,  "[.]", 4)
    genes2 <- ._split_and_pad(g_second_s, "[.]", 4)
    
    # Remap to V/D/J/C (D=NA for 3-gene chains)
    g1 <- .remap_VDJC(genes1)
    g2 <- .remap_VDJC(genes2)
    
    # Display helpers
    pair_strict <- function(lhs, rhs) {
      lhs <- ifelse(is.na(lhs) | lhs == "", "NA", as.character(lhs))
      rhs <- ifelse(is.na(rhs) | rhs == "", "NA", as.character(rhs))
      paste(lhs, rhs, sep = ";")
    }
    
    # For genes: if chain absent -> "NA"
    show_gene <- function(g, present) {
      out <- ifelse(is.na(g) | g == "", "NA", g)
      out[!present] <- "NA"
      out
    }
    # For D gene: if chain present but D missing (3-gene) -> "None"; if chain absent -> "NA"
    show_D <- function(d, present) {
      out <- d
      out <- ifelse(is.na(out) | out == "", "None", out)  
      out[!present] <- "NA"                               
      out
    }
    
    # Build display vectors
    V1 <- show_gene(g1$V, present1); V2 <- show_gene(g2$V, present2)
    D1 <- show_D(  g1$D, present1); D2 <- show_D(  g2$D, present2)
    J1 <- show_gene(g1$J, present1); J2 <- show_gene(g2$J, present2)
    C1 <- show_gene(g1$C, present1); C2 <- show_gene(g2$C, present2)
    
    # CDR3: if chain absent or blank -> "NA"
    aa1 <- ifelse(present1 & !(is.na(aa_first)  | aa_first  == ""), aa_first,  "NA")
    aa2 <- ifelse(present2 & !(is.na(aa_second) | aa_second == ""), aa_second, "NA")
    nt1 <- ifelse(present1 & !(is.na(nt_first)  | nt_first  == ""), nt_first,  "NA")
    nt2 <- ifelse(present2 & !(is.na(nt_second) | nt_second == ""), nt_second, "NA")
    
    # Assemble; always "lhs;rhs"
    out <- data.frame(
      Clones      = result[["Clones"]],
      Proportion  = result[["Proportion"]],
      CDR3.nt     = pair_strict(nt1, nt2),
      CDR3.aa     = pair_strict(aa1, aa2),
      V.name      = pair_strict(V1, V2),
      D.name      = pair_strict(D1, D2),  # preserves "None" vs "NA"
      J.name      = pair_strict(J1, J2),
      C.name      = pair_strict(C1, C2),
      Barcode     = result[["barcode"]],
      stringsAsFactors = FALSE
    )
    out
  })
  
  names(data_out) <- names(df_list)
  list(data = data_out, meta = meta)
}