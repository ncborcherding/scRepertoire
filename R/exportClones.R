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
#' @param df The product of \code{\link{combineTCR}}, \code{\link{combineBCR}}, or
#'  \code{\link{combineExpression}}.
#' @param group.by The column header used for grouping.
#' @param split.by If using a single-cell object, the column header 
#' to group the new list. NULL will return clusters.
#' @param write.file TRUE, save the file, FALSE, return a data.frame
#' @param dir directory location to save the csv
#' @param file.name the csv file name
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split
#' @importFrom utils write.csv
#' @export
#' @return ggplot of percentage of V and J gene pairings as a heatmap
#' @author Jonathan Noonan, Nick Borcherding
exportClones <- function(df,
                         group.by = NULL,
                         split.by = NULL,
                         write.file = TRUE,
                         dir = NULL, 
                         file.name = "clones.csv") {
  
  df <- .data.wrangle(df, split.by, "CTgene", "both")

  df <- bind_rows(df, .id = "group")
  
  genes <- str_split(df[,"CTgene"], "_", simplify = TRUE)
  aa <- str_split(df[,"CTaa"], "_", simplify = TRUE)
  nt <- str_split(df[,"CTnt"], "_", simplify = TRUE)
  chain1_gene <- str_split(genes[,1], "[.]", simplify = TRUE)
  
  mat <- data.frame(row.names = df[,"barcode"], 
                    chain1_aa = aa[,1], 
                    chain1_nt = nt[,1], 
                    chain1_genes = genes[,1], 
                    chain2_aa = aa[,2],
                    chain2_nt = nt[,2], 
                    chain2_genes = genes[,2],
                    group = df[,"group"])
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