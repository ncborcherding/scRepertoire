#' Deconvolute Contig Information from Multiplexed Experiments
#'
#' This function reprocess and forms a list of contigs for downstream analysis 
#' in scRepertoire, [createHTOContigList()] take the filtered contig 
#' annotation output and the single-cell RNA object to create the list. 
#' If using an integrated single-cell object, it is recommended to split the 
#' object by sequencing run and remove extra prefixes and suffixes on the 
#' barcode before using [createHTOContigList()]. Alternatively, 
#' the variable **multi.run** can be used to separate a list of contigs
#' by a meta data variable. This may have issues with the repeated barcodes.
#' 
#' @examples
#' \dontrun{
#' filtered.contig <- read.csv(".../Sample/outs/filtered_contig_annotations.csv")
#' 
#' contig.list <- createHTOContigList(contig = filtered.contig, 
#'                                    sc.data = Seurat.Obj, 
#'                                    group.by = "HTO_maxID")
#' }
#' 
#' @param contig The filtered contig annotation file from multiplexed experiment
#' @param sc.data The Seurat or Single-Cell Experiment object.
#' @param group.by One or more meta data headers to create the contig 
#' list based on. If more than one header listed, the function combines 
#' them into a single variable.
#' @param multi.run If using integrated single-cell object, the meta data 
#' variable that indicates the sequencing run.
#' @export
#' @concept Loading_and_Processing_Contigs
#' @return Returns a list of contigs as input for [combineBCR()] 
#' or [combineTCR()]


createHTOContigList <- function(contig, 
                                sc.data, 
                                group.by = NULL, 
                                multi.run = NULL){
  contig.list <- NULL
  .checkSingleObject(sc.data)
  meta <- .grabMeta(sc.data)
  if (length(group.by) > 1) {
    meta["group.by"] <- apply(meta[ , group.by] , 1 , paste , collapse = "." )
  } else {
    meta[,"group.by"] <- meta[ , group.by]
  }
  unique.groups <- unique(meta$group.by)
  
  if (!is.null(multi.run) & inherits(x=contig, what ="list")) {
    tmp.list <- list()
    runs <- unique(meta[,multi.run])
    for (x in seq_along(contig)) {
      cont.tmp <- contig[[x]]
      cont.tmp$barcode <- paste0(cont.tmp$barcode, "_", x)
      meta.tmp <- meta[meta[,multi.run] == runs[x],] #Subset meta in batch
      rownames(meta.tmp) <- substring(rownames(meta.tmp), nchar(rownames(meta.tmp)[1])-17, nchar(rownames(meta.tmp)[1]))
      cont.tmp <- cont.tmp[cont.tmp$barcode %in% row.names(meta.tmp) ] #subset contig in new meta
      for (y in seq_along(unique.groups)) {
        sub.con <- cont.tmp[cont.tmp$barcode %in% rownames(subset(meta.tmp, group.by == unique.groups[y])),]
        tmp.list[[y]] <- sub.con
      }
      if(x == 1) {
        contig.list  <- tmp.list
      } else {
        contig.list <- c(contig.list, tmp.list)
      }
    }
  } else {
    cont.tmp <- contig[contig$barcode %in% rownames(meta), ]
    for (i in seq_along(unique.groups)) {
      sub.con <- cont.tmp[cont.tmp$barcode %in% rownames(subset(meta, group.by == unique.groups[i])),]
      contig.list[[i]] <- sub.con
    }
    names(contig.list) <- unique.groups
  }
  contig.list
}
