#' Exporting clonal information
#'
#' This function saves a csv file of clones (genes, amino acid, and 
#' nucleotide sequences) by barcodes. **format** determines 
#' the structure of the csv file.
#' \itemize{
#'   \item{"paired"} will export sequences by barcodes and include multiple chains.
#'   \item{"airr"} will export a data frame that is consistent with the AIRR format.
#'   \item{"TCRMatch"} will export a data frame that has the TRB chain with count information.
#'   \item{"tcrpheno"} will export a data frame with TCR-alpha and TCR-beta chains in separate columns.
#' }
#' 
#' @examples
#' \dontrun{
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' exportClones(combined, 
#'              format = "paired")
#' }
#'                                    
#' @param input.data The product of [combineTCR()], 
#' [combineBCR()], or [combineExpression()].
#' @param format The format to export the clones - "paired", "airr", "TCRMatch",
#' or "tcrpheno".
#' @param group.by The variable to use for grouping.
#' @param write.file **TRUE**, save the file or **FALSE**, 
#' return a data.frame
#' @param dir directory location to save the csv
#' @param file.name the csv file name
#' @importFrom utils write.csv
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return CSV file of the paired sequences.
#' @author Jonathan Noonan, Nick Borcherding
exportClones <- function(input.data,
                         format = "paired",
                         group.by = NULL,
                         write.file = TRUE,
                         dir = NULL, 
                         file.name = "clones.csv") {
  
  summaryFunc <- switch(format,
                       "paired"     = .PairedExport,
                       "airr"       = .AIRRExport,
                       "TCRMatch"   = .TCRmatchExport,
                       "tcrpheno"   = .tcrphenoExport,
                       stop("Invalid format provided. Options are: 'paired', 'airr', 'TCRMatch', 'tcrpheno'"))
  mat <- summaryFunc(input.data)
  
  mat[mat == "NA"] <- NA
  if(!write.file) {
    return(mat)
  }
  
  if(is.null(dir)) {
    dir <- "."
  }
  filepath <- file.path(dir, file.name)
  write.csv(mat, file = filepath)
}


#' @importFrom stats na.omit
.TCRmatchExport <- function(input.data, group.by) {
  input.data <- .data.wrangle(input.data, NULL, "CTgene", "TRB")
  
  for(i in seq_along(input.data)) {
    input.data[[i]] <- .off.the.chain(input.data[[i]], "TRB", "CTaa", check = FALSE)
    input.data[[i]] <- .off.the.chain(input.data[[i]], "TRB", "CTnt", check = FALSE)
  }
  
  list_names <- names(input.data)
  if (is.null(list_names)) {
    list_names <- as.character(seq_along(input.data))
  }
  input.data <- mapply(function(df, name) transform(df, group = name), 
                       input.data, list_names, SIMPLIFY = FALSE)
  input.data <- do.call(rbind, input.data)
  
  #Remove cells without TRB
  if(any(c("IGH", "IGL", "IGK", "TRG", "TRA", "TRD") %in% unique(substr(input.data$CTgene, 1,3)))) {
    input.data <-  input.data[-grep("IGH|TRG|IGL|IGK|TRA|TRD", substr(input.data$CTgene, 1,3)),]
  }
  
  mat <- data.frame(chain2_aa = substring(input.data[,"CTaa"], 2, nchar(input.data[,"CTaa"])),
                    group = input.data[,"group"])
  mat <- na.omit(mat)
  
  mat$clonalFrequency <- ave(rep(1, nrow(mat)), mat[,c("group", "chain2_aa")], FUN = sum)
  mat <- unique(mat)
  
  return(mat)
}

.PairedExport <- function(input.data, group.by) {
  input.data <- .data.wrangle(input.data, group.by, "CTgene", "both")
  
  list_names <- names(input.data)
  if (is.null(list_names)) {
    list_names <- as.character(seq_along(input.data))
  }
  input.data <- mapply(function(df, name) transform(df, group = name), 
                       input.data, list_names, SIMPLIFY = FALSE)
  input.data <- do.call(rbind, input.data)
  
  str_to_matrix <- function(x) {
    s <- strsplit(x, "_")
    t(sapply(s, `[`, 1:2))
  }
  genes <- str_to_matrix(input.data[,"CTgene"])
  aa <- str_to_matrix(input.data[,"CTaa"])
  nt <- str_to_matrix(input.data[,"CTnt"])
  
  mat <- data.frame(row.names = input.data[,"barcode"], 
                    chain1_aa = aa[,1], 
                    chain1_nt = nt[,1], 
                    chain1_genes = genes[,1], 
                    chain2_aa = aa[,2],
                    chain2_nt = nt[,2], 
                    chain2_genes = genes[,2],
                    group = input.data[,"group"])
  return(mat)
}

.AIRRExport <- function(input.data, group.by) {
  input.data <- .data.wrangle(input.data, group.by, "CTgene", "both")
  
  list_names <- names(input.data)
  if (is.null(list_names)) {
    list_names <- as.character(seq_along(input.data))
  }
  input.data <- mapply(function(df, name) transform(df, group = name), 
                       input.data, list_names, SIMPLIFY = FALSE)
  input.data <- do.call(rbind, input.data)
  
  mat <- list()
  
  str_to_matrix <- function(x, split) {
    s <- strsplit(x, split)
    t(sapply(s, `[`, 1:max(sapply(s, length))))
  }
  
  .process_row <- function(row) {
    genes <- str_to_matrix(row$CTgene, "_")
    aa <- str_to_matrix(row$CTaa, "_")
    nt <- str_to_matrix(row$CTnt, "_")
    
    if (any(grepl(";", aa))) {
      multi.chain.pos <- grep(";", aa)
      genes[, multi.chain.pos] <- sapply(genes[, multi.chain.pos], function(x) strsplit(x, ";")[[1]][1])
      aa[, multi.chain.pos] <- sapply(aa[, multi.chain.pos], function(x) strsplit(x, ";")[[1]][1])
      nt[, multi.chain.pos] <- sapply(nt[, multi.chain.pos], function(x) strsplit(x, ";")[[1]][1])
    }
    
    chain1_gene <- str_to_matrix(genes[, 1], "[.]")
    chain2_gene <- str_to_matrix(genes[, 2], "[.]")
    locus1 <- substr(chain1_gene[, 1], 1, 3)
    locus2 <- substr(chain2_gene[, 1], 1, 3)
    
    .sort_gene_calls <- function(gene) {
      if (ncol(gene) == 3) {
        gene <- cbind(gene, NA)
        colnames(gene) <- c("v_call", "j_call", "c_call", "d_call")
      } else if (ncol(gene) == 4) {
        colnames(gene) <- c("v_call", "d_call", "j_call", "c_call")
      } else if (all(is.na(gene))) {
        gene <- matrix(ncol = 4, nrow = 1, NA)
        colnames(gene) <- c("v_call", "j_call", "c_call", "d_call")
      }
      return(gene)
    }
    
    chain1_gene <- .sort_gene_calls(chain1_gene)
    chain2_gene <- .sort_gene_calls(chain2_gene)
    
    tmp.out <- data.frame(
      cell_id = row$barcode,
      locus = c(locus1, locus2),
      v_call = c(chain1_gene[,"v_call"], chain2_gene[,"v_call"]),
      d_call = c(chain1_gene[,"d_call"], chain2_gene[,"d_call"]),
      j_call = c(chain1_gene[,"j_call"], chain2_gene[,"j_call"]),
      c_call = c(chain1_gene[,"c_call"], chain2_gene[,"c_call"]),
      junction = c(nt[, 1], nt[, 2]),
      junction_aa = c(aa[, 1], aa[, 2])
    )
    
    tmp.out[tmp.out == "NA"] <- NA
    na.to.remove <- which(rowSums(is.na(tmp.out)) > 4)
    if (length(na.to.remove) >= 1) {
      tmp.out <- tmp.out[-na.to.remove, ]
    }
    
    return(tmp.out)
  }
  
  mat <- lapply(seq_len(nrow(input.data)), function(i) .process_row(input.data[i,]))
  mat <- do.call(rbind, mat)
  return(mat)
}

.tcrphenoExport <- function(input.data, group.by) {
  input.data <- .data.wrangle(input.data, NULL, "CTgene", "both")
  dat <- do.call(rbind, unname(input.data))
  
  # split strings into a matrix with a fixed number of columns
  split_to_fixed_matrix <- function(x, split_char, num_cols) {
    s <- strsplit(x, split_char)
    t(sapply(s, `[`, seq_len(num_cols)))
  }
  
  aa_split <- split_to_fixed_matrix(dat$CTaa, "_", 2)
  nt_split <- split_to_fixed_matrix(dat$CTnt, "_", 2)
  gene_split <- split_to_fixed_matrix(dat$CTgene, "_", 2)
  
  tcr1_genes <- split_to_fixed_matrix(gene_split[, 1], "[.]", 2)
  tcr2_genes <- split_to_fixed_matrix(gene_split[, 2], "[.]", 3)[,c(1,3)]
  
  contigs <- data.frame(
    cell = rownames(dat),
    TCRA_cdr3aa = aa_split[, 1],
    TCRA_vgene = tcr1_genes[, 1],
    TCRA_jgene = tcr1_genes[, 2],
    TCRA_cdr3nt = nt_split[, 1],
    TCRB_cdr3aa = aa_split[, 2],
    TCRB_vgene = tcr2_genes[, 1],
    TCRB_jgene = tcr2_genes[, 2],
    TCRB_cdr3nt = nt_split[, 2]
  )
  # Clean multichains
  columns_to_clean <- c("TCRA_cdr3aa", "TCRA_cdr3nt", "TCRB_cdr3aa", "TCRB_cdr3nt")
  contigs[columns_to_clean] <- lapply(contigs[columns_to_clean], function(x) {
    gsub(";.*", "", x)
  })
  
  contigs[contigs == ""] <- NA
  contigs[contigs == "NA"] <- NA
  
  return(contigs)
}
