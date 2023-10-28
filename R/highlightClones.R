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
#' scRep_example  <- combineExpression(combined, 
#'                                     scRep_example)
#' 
#' #Using highlightClonotype()
#' scRep_example   <- highlightClones(scRep_example, 
#'                                    cloneCall= "aa", 
#'                                    sequence = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF"))
#' 
#' @param sc.data The single-cell object to attach after \code{\link{combineExpression}}
#' @param cloneCall How to call the clonotype - VDJC gene (gene), 
#' CDR3 nucleotide (nt), CDR3 amino acid (aa),
#' VDJC gene + CDR3 nucleotide (strict) or a custom variable in the data. 
#' @param sequence The specific sequence or sequence to highlight
#' @importFrom S4Vectors DataFrame
#' @export
#' @concept SC_Functions
#' @return Single-cell object object with new meta data column for indicated clones
highlightClones <- function(sc.data, 
                            cloneCall = c("gene", "nt", "aa", "strict"), 
                            sequence = NULL){
  if (!is_seurat_or_se_object(sc.data)) {
    stop("Please select a single-cell object") 
  }
  
  cloneCall <- .theCall(.grabMeta(sc.data), cloneCall)
  meta <- .grabMeta(sc.data)
  meta$highlight <- NA
  for(i in seq_along(sequence)) {
    meta$highlight <-  ifelse(meta[,cloneCall] == sequence[i], 
                              sequence[i], meta$highlight)
  }
  
  meta <- meta[,-(which(colnames(meta) == "ident"))]
  
  if(is_se_object(sc.data)) {
    colData(sc.data) <- DataFrame(meta)
  } else {
    col.name <- names(meta) %||% colnames(meta)
    sc.data[[col.name]] <- meta
  }
  return(sc.data)
}