#' Highlighting Specific Clones
#'
#' Use a specific clonal sequence to highlight on top of the dimensional 
#' reduction in single-cell object.
#'
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example  <- get(data("scRep_example"))
#' 
#' # Using combineExpresion()
#' scRep_example  <- combineExpression(combined, 
#'                                     scRep_example)
#' 
#' # Using highlightClones()
#' scRep_example   <- highlightClones(scRep_example, 
#'                                    cloneCall= "aa", 
#'                                    sequence = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF"))
#' 
#' @param sc.data The single-cell object to attach after 
#' [combineExpression()]
#' @param cloneCall Defines the clonal sequence grouping. Accepted values 
#' are: `gene` (VDJC genes), `nt` (CDR3 nucleotide sequence), `aa` (CDR3 amino 
#' acid sequence), or `strict` (VDJC). A custom column header can also be used.
#' @param sequence The specific sequence or sequence to highlight
#' @importFrom S4Vectors DataFrame
#' @export
#' @concept SC_Functions
#' @return Single-cell object object with new meta data column 
#' for indicated clones
highlightClones <- function(sc.data, 
                            cloneCall = c("gene", "nt", "aa", "strict"), 
                            sequence = NULL){
  if (!.is.seurat.or.se.object(sc.data)) {
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
  
  if(.is.se.object(sc.data)) {
    colData(sc.data) <- DataFrame(meta)
  } else {
    col.name <- names(meta) %||% colnames(meta)
    sc.data[[col.name]] <- meta
  }
  return(sc.data)
}
