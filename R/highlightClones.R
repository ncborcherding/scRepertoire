#' Highlighting specific clonotypes in Seurat
#'
#' Use a specific clonotype sequence to highlight on top of the dimensional 
#' reduction in single-cell object.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' 
#' #Using combineExpresion()
#' scRep_example  <- combineExpression(combined, scRep_example)
#' 
#' #Using highlightClonotype()
#' scRep_example   <- highlightClones(scRep_example , cloneCall= "aa", 
#' sequence = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF"))
#' 
#' @param sc The single-cell object to attach after \code{\link{combineExpression}}
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa), or 
#' VDJC gene + CDR3 nucleotide (strict).
#' @param sequence The specific sequence or sequence to highlight
#' @importFrom S4Vectors DataFrame
#' @export
#' @return Single-cell object object with new meta data column for indicated clones
highlightClones <- function(sc, 
                            cloneCall = c("gene", "nt", "aa", "strict"), 
                            sequence = NULL){
  if (!is_seurat_or_se_object(sc)) {
    stop("Please select a single-cell object") 
  }
  
  cloneCall <- .theCall(cloneCall)
  meta <- .grabMeta(sc)
  meta$highlight <- NA
  for(i in seq_along(sequence)) {
    meta$highlight <-  ifelse(meta[,cloneCall] == sequence[i], 
                              sequence[i], meta$highlight)
  }
  
  meta <- meta[,-(which(colnames(meta) == "ident"))]
  
  if(is_se_object(sc)) {
    colData(sc) <- DataFrame(meta)
  } else {
    col.name <- names(meta) %||% colnames(meta)
    sc[[col.name]] <- meta
  }
  return(sc)
}