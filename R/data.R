#' A data set of T cell contigs as a list outputed from the 
#' filter_contig_annotation files.
#' @docType data
#' @name contig_list
#' 
NULL

#' A seurat object of 500 single T cells derived
#' from 4 patients with acute respiratory distress (ARDS).
#' 
#' @description The object is compatible with `contig_list` and the TCR
#' sequencing data can be added with `combineExpression`.
#' 
#' @name scRep_example
#' @docType data
#'
NULL

#' Processed subset of `contig_list`
#' 
#' @description A list of 8 dataframes of T cell contigs outputted from the
#' `filtered_contig_annotation` files, but subsetted to  365 valid T cells
#' which correspond to the same barcodes found in `scRep_example`
#'
#' @usage data("mini_contig_list")
#'
#' @format An R `list` of `data.frame` objects
#' 
#' @docType data
#'
#' @seealso \code{\link{contig_list}}
#'
"mini_contig_list"
