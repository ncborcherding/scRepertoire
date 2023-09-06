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
#' @name scRep_example
#' @docType data
#'
NULL

#' Processed subset of `contig_list`
#' 
#' @description A list of 6 dataframes of T cell contigs outputted from the
#' `filtered_contig_annotation` files, but subsetted to about 92 valid T cells
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

# # Code used for creating the mini_contig_list:

# library(hash, usethis)
# 
#data("contig_list", "scRep_example")

#mini_contig_list <- combineTCR(
#	contig_list,
#	samples = c("PY", "PY", "PX", "PX", "PZ","PZ"),
#	ID = c("P", "T", "P", "T", "P", "T")
#)
#all_barcodes <- names(scRep_example@active.ident)
#barcode_set <- hash::hash(all_barcodes, all_barcodes) # a worse version of a set
#col_names <- names(mini_contig_list[[1]])

#for (i in seq_along(mini_contig_list)) {
#	curr_df <- setNames(
#		data.frame(replicate(length(col_names), character(0))), col_names
#	)
#	len <- 0
#	for (j in seq_along(mini_contig_list[[i]][[1]])) {
#		if (is.null(barcode_set[[mini_contig_list[[i]][[1]][[j]]]])) {
#			next
#		}
#		len <- len + 1
#		curr_df[len, ] <- mini_contig_list[[i]][j, ]
#	}
#	mini_contig_list[[i]] <- curr_df
#}
#usethis::use_data(mini_contig_list)
