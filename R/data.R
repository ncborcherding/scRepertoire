#' A data set of T cell contigs as a list outputed from the 
#' filter_contig_annotation files.
#' @docType data
#' @name contig_list
#' 
NULL

#' A seurat object of 100 single T cells derived
#' from 3 clear cell renal carcinoma patients.
#' 
#' @description The object is compatible with `contig_list` and the TCR
#' sequencing data can be added with `combineExpression`.
#' 
#' @name screp_example
#' @docType data
#'
NULL

#' Processed subset of `contig_list`
#' 
#' @description A list of 6 dataframes of T cell contigs outputted from the
#' `filtered_contig_annotation` files, but subsetted to about 92 valid T cells
#' which correspond to the same barcodes found in `screp_example`
#'
#' @usage data("combined_mini_contig_list")
#'
#' @format An R `list` of `data.frame` objects
#' 
#' @docType data
#'
#' @seealso \code{\link{contig_list}}
#'
"combined_mini_contig_list"
