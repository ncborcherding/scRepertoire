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

# # Code used for creating the combined_mini_contig_list:

# library(hash, usethis)
# 
#data("contig_list", "screp_example")

#combined_mini_contig_list <- combineTCR(
#	contig_list,
#	samples = c("PY", "PY", "PX", "PX", "PZ","PZ"),
#	ID = c("P", "T", "P", "T", "P", "T")
#)
#all_barcodes <- names(screp_example@active.ident)
#barcode_set <- hash::hash(all_barcodes, all_barcodes) # a worse version of a set
#col_names <- names(combined_mini_contig_list[[1]])

#for (i in seq_along(combined_mini_contig_list)) {
#	curr_df <- setNames(
#		data.frame(replicate(length(col_names), character(0))), col_names
#	)
#	len <- 0
#	for (j in seq_along(combined_mini_contig_list[[i]][[1]])) {
#		if (is.null(barcode_set[[combined_mini_contig_list[[i]][[1]][[j]]]])) {
#			next
#		}
#		len <- len + 1
#		curr_df[len, ] <- combined_mini_contig_list[[i]][j, ]
#	}
#	combined_mini_contig_list[[i]] <- curr_df
#}
#usethis::use_data(combined_mini_contig_list)
