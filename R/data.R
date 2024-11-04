#' A list of 8 single-cell T cell receptor sequences runs.
#' 
#' @description A list of 8 `filtered_contig_annotations.csv` files
#' outputted from 10X Cell Ranger. More information on the 
#' data can be found in the following
#'  [manuscript](https://pubmed.ncbi.nlm.nih.gov/33622974/).
#'  
#' @docType data
#' @concept Data
#' @name contig_list
#' 
NULL

#' A Seurat object of 500 single T cells,
#' 
#' @description The object is compatible with `contig_list` and the TCR
#' sequencing data can be added with `combineExpression`.  The data is 
#' from 4 patients with acute respiratory distress, with samples taken
#' from both the lung and peripheral blood. More information on the 
#' data can be found in the following
#'  [manuscript](https://pubmed.ncbi.nlm.nih.gov/33622974/).
#' 
#' @name scRep_example
#' @concept Data
#' @docType data
#'
NULL

#' Processed subset of `contig_list`
#' 
#' @description A list of 8 data frames of T cell contigs outputted from the
#' `filtered_contig_annotation` files, but subsetted to  365 valid T cells
#' which correspond to the same barcodes found in `scRep_example`. The
#' data is originally derived from the following
#'  [manuscript](https://pubmed.ncbi.nlm.nih.gov/33622974/).
#'
#' @usage data("mini_contig_list")
#'
#' @format An R `list` of `data.frame` objects
#' 
#' @concept Data
#' 
#' @docType data
#'
#' @seealso [contig_list()]
#'
"mini_contig_list"
