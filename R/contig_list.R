#' A data set T cell contigs as a list outputed from the filter_contig_annotation files.
#' @name contig_list
#' @docType list
#'
#' \itemize{
#' \item barcode. the cell barcodes for sequenced cells
#' \item is_cell. A binary classifier
#' \item contig_id. Cellranger contig assignment
#' \item high_confidence. Binary classifier of if the chain sequence is assigned with high confidence
#' \item length.
#' \item v_gene. corresponding V gene of the immune receptor
#' \item d_gene. corresponding D gene of the immune receptor
#' \item j_gene. corresponding J gene of the immune receptor
#' \item C_gene. corresponding C gene of the immune receptor
#' \item full_length.
#' \item productive. Binary classifier if the chain is estimated to be expressed
#' \item cdr3. The amino acid sequence of the chain
#' \item cd3_nt. The nucleotide sequence of the read
#' \item reads. The count of sequences assigned to the read
#' \item umis. the unique molecular identifiers
#' \item raw_clonotype_id. Cellranger colonotype assignment
#' \item raw_consensus_id. Cellranger consensus clonotype assignment
#' }
#'
#' @source https://github.com/ncborcherding/scRepertoire/blob/master/data/contig_list.rda
NULL
