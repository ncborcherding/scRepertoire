#' Exporting clonotypes
#'
#' This function saves a csv file of clonotypes (genes, amino acid, and nucleotide sequences)
#' by barcodes.
#' 
#' @examples
#' \dontrun{
#' #Making combined contig data
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' exportClones(combined)
#' }
#'                                    
#' @param input.data The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param group.by The variable to use for grouping.
#' @param write.file \strong{TRUE}, save the file or \strong{FALSE}, return a data.frame
#' @param dir directory location to save the csv
#' @param file.name the csv file name
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split
#' @importFrom utils write.csv
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return CSV file of the paired sequences.
#' @author Jonathan Noonan, Nick Borcherding
exportClones <- function(input.data,
                         group.by = NULL,
                         write.file = TRUE,
                         dir = NULL, 
                         file.name = "clones.csv") {
  
  input.data <- .data.wrangle(input.data, group.by, "CTgene", "both")

  input.data <- bind_rows(input.data, .id = "group")
  
  genes <- str_split(input.data[,"CTgene"], "_", simplify = TRUE)
  aa <- str_split(input.data[,"CTaa"], "_", simplify = TRUE)
  nt <- str_split(input.data[,"CTnt"], "_", simplify = TRUE)
  chain1_gene <- str_split(genes[,1], "[.]", simplify = TRUE)
  
  mat <- data.frame(row.names = input.data[,"barcode"], 
                    chain1_aa = aa[,1], 
                    chain1_nt = nt[,1], 
                    chain1_genes = genes[,1], 
                    chain2_aa = aa[,2],
                    chain2_nt = nt[,2], 
                    chain2_genes = genes[,2],
                    group = input.data[,"group"])
  mat[mat == "NA"] <- NA
  if(!write.file) {
    return(mat)
  }
  
  if(is.null(dir)) {
    dir <- "./"
  }
  filepath <- paste0(dir, file.name)
  write.csv(mat, file = filepath)
}