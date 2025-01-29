#' Exporting clones
#'
#' This function saves a csv file of clones (genes, amino acid, and 
#' nucleotide sequences) by barcodes. **format** determines 
#' the structure of the csv file - *paired* will export sequences 
#' by barcodes and include multiple chains, *airr* will export a data 
#' frame that is consistent with the AIRR format, and *TCRMatch* will 
#' export a data frame that has the TRB chain with count information.
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
#' @param format The format to export the clones - "paired", "airr", or "TCRMatch".
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
                      "paired" = .PairedExport,
                      "airr"  = .AIRRExport,
                      "TCRMatch" = .TCRmatchExport,
                      stop("Invalid format provided"))
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


.TCRmatchExport <- function(input.data) {

  input.data <- .data.wrangle(input.data, NULL, "CTgene", "TRB")
  
  for(i in seq_along(input.data)) {
    input.data[[i]] <- .off.the.chain(input.data[[i]], "TRB", "CTaa", check = FALSE)
    input.data[[i]] <- .off.the.chain(input.data[[i]], "TRB", "CTnt", check = FALSE)
  }
  
  input.data <- bind_rows(input.data, .id = "group")
  #Remove cells without TRB
  if(any(c("IGH", "IGL", "IGK", "TRG", "TRA", "TRD") %in% unique(substr(input.data$CTgene, 1,3)))) {
    input.data <-  input.data[-grep("IGH|TRG|IGL|IGK|TRA|TRD", substr(input.data$CTgene, 1,3)),]
  }
  
  mat <- data.frame(chain2_aa = substring(input.data[,"CTaa"], 2, nchar(input.data[,"CTaa"])),
                    group = input.data[,"group"])
  mat <- na.omit(mat)
  mat <- mat %>%
          group_by(group, chain2_aa) %>%
          dplyr::mutate(clonalFrequency = n())
  mat <- unique(mat)
  
  return(mat)
}

.PairedExport <- function(input.data) {
  input.data <- .data.wrangle(input.data, group.by, "CTgene", "both")
  
  input.data <- bind_rows(input.data, .id = "group")
  
  genes <- str_split(input.data[,"CTgene"], "_", simplify = TRUE)
  aa <- str_split(input.data[,"CTaa"], "_", simplify = TRUE)
  nt <- str_split(input.data[,"CTnt"], "_", simplify = TRUE)
  
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

.AIRRExport <- function(input.data) {
  input.data <- .data.wrangle(input.data, group.by, "CTgene", "both")
  input.data <- bind_rows(input.data, .id = "group")
  
  mat <- list()

    
    mat <- list()
    
    # Define a function for processing each row - main speed bottleneck
    .process_row <- function(row) {
      # Split gene, amino acid, and nucleotide columns
      genes <- str_split(row$CTgene, "_", simplify = TRUE)
      aa <- str_split(row$CTaa, "_", simplify = TRUE)
      nt <- str_split(row$CTnt, "_", simplify = TRUE)
      
      # Remove secondary chain if multiple chains
      if (any(grepl(";", aa))) {
        multi.chain.pos <- grep(";", aa)
        genes[, multi.chain.pos] <- sapply(genes[, multi.chain.pos], function(x) 
          str_split(x, ";", simplify = TRUE)[, 1])
        aa[, multi.chain.pos] <- sapply(aa[, multi.chain.pos], function(x) 
          str_split(x, ";", simplify = TRUE)[, 1])
        nt[, multi.chain.pos] <- sapply(nt[, multi.chain.pos], function(x) 
          str_split(x, ";", simplify = TRUE)[, 1])
      }
      
      # Extract locus information
      chain1_gene <- str_split(genes[, 1], "[.]", simplify = TRUE)
      chain2_gene <- str_split(genes[, 2], "[.]", simplify = TRUE)
      locus1 <- substr(chain1_gene[, 1], 1, 3)
      locus2 <- substr(chain2_gene[, 1], 1, 3)
      
      # Define a function to sort gene calls
      .sort_gene_calls <- function(gene) {
        if (length(gene) == 3) {
          gene <- cbind(gene, NA)
          colnames(gene) <- c("v_call", "j_call", "c_call", "d_call")
        } else if (length(gene) == 4) {
          colnames(gene) <- c("v_call", "d_call", "j_call", "c_call")
        } else if (all(gene == "NA")) {
          gene <- matrix(ncol = 4, nrow = 1, NA)
          colnames(gene) <- c("v_call", "j_call", "c_call", "d_call")
        }
        return(gene)
      }
      
      # Sort gene calls for both chains
      chain1_gene <- .sort_gene_calls(chain1_gene)
      chain2_gene <- .sort_gene_calls(chain2_gene)
      
      # Create the formatted data frame - smaller speed bottleneck of ~10 sec
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
      
      # Replace "NA" with NA in the data frame
      tmp.out[tmp.out == "NA"] <- NA
      
      # Identify rows with more than 4 NA values and remove them
      na.to.remove <- which(rowSums(is.na(tmp.out)) > 4)
      if (length(na.to.remove) >= 1) {
        tmp.out <- tmp.out[-na.to.remove, ]
      }
      
      return(tmp.out)
    }
    
    # Process each row and store the results in the list
    for (i in seq_len(nrow(input.data))) {
      mat[[i]] <- .process_row(input.data[i, ]) # main speed bottleneck (see above)
    }
  mat <- do.call(rbind, mat)
  return(mat)
}
