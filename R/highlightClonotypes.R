#' Highlighting specific clonotypes in Seurat
#'
#' Use a specific clonotype sequence to highlight on top of the dimensional 
#' reduction in seurat object.
#'
#' @examples
#' #' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example "))
#' 
#' #Using combineExpresion()
#' scRep_example  <- combineExpression(combined, scRep_example  )
#' 
#' #Using highlightClonotype()
#' scRep_example   <- highlightClonotypes(scRep_example , cloneCall= "aa", 
#' sequence = c("CAVNGGSQGNLIF_CSAEREDTDTQYF"))
#' 
#' @param sc The Seurat object to attach
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param sequence The specific sequence or sequence to highlight
#'
#' @export
#' @return Seurat object with new meta data column for indicated clones
highlightClonotypes <- function(sc, 
                                cloneCall = c("gene", "nt", "aa", "strict"), 
                                sequence = NULL){
  if (!inherits(x=sc, what ="Seurat")) {
    stop("Object indicated is not of class 'Seurat', make sure you 
            are using the correct data.") }
  cloneCall <- theCall(cloneCall)
  meta <- grabMeta(sc)
  meta$highlight <- NA
  for(i in seq_along(sequence)) {
    meta$highlight <-  ifelse(meta[,cloneCall] == sequence[i], 
                              sequence[i], meta$highlight)
  }
  meta <- meta[,-(which(colnames(meta) == "ident"))]
  col.name <- names(meta) %||% colnames(meta)
  sc[[col.name]] <- meta
  sc
}