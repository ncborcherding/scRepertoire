#' Annotate invariant T cells (MAIT or iNKT) in single-cell TCR data
#' 
#' The [scoreInvariant()] function identifies mucosal-associated invariant T (MAIT) 
#' cells or invariant natural killer T (iNKT) cells from single-cell sequencing 
#' datasets based on their characteristic T-cell receptor (TCR) usage. It 
#' extracts TCR chain information from the provided single-cell data, checks 
#' it against known invariant T-cell receptor criteria for either MAIT or iNKT 
#' cells, and returns a score indicating the presence (1) or absence (0) of 
#' these invariant cell populations for each individual cell. The function 
#' supports data from mouse and human samples, providing a convenient method 
#' to annotate specialized T-cell subsets within single-cell analyses.
#'
#' @param input.data The single-cell object after [combineExpression()]
#' @param type Character specifying the type of invariant cells to
#'annotate ('MAIT' or 'iNKT').
#' @param species Character specifying the species 
#' ('mouse' or 'human').
#' 
#' @return A data.frame with barcode identifiers and corresponding annotation 
#' scores (0 or 1).
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' #Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' #Using annotateInvariant
#' annotateInvariant(input.data = seurat.obj, type = "MAIT", species = "human")
#' annotateInvariant(input.data = seurat.obj, type = "iNKT", species = "human")
#' 
#' @importFrom immApex getIR
#'
#' @export
annotateInvariant <- function(input.data, 
                              type = c("MAIT", "iNKT"), 
                              species = c("mouse", "human")) {
  
  if(!is_seurat_or_se_object(input.data)) {
    stop("Please ensure input.data is a single-cell object after 
         combineExpression()")
  }
  
  type <- match.arg(type)
  species <- match.arg(species)
  
  TCRS <- lapply(c("TRA", "TRB") function(x) {
    tmp <- getIR(input.data, chains = x, sequence.type = "aa")
    tmp
  })
  
  criteria <- switch(type,
                     "MAIT" = .MAIT.criteria,
                     "iNKT" = .iNKT.criteria)
  
  species.criteria <- criteria[[species]]
  
  TRA.data <- TCRS[[1]]
  TRB.data <- TCRS[[2]]
  
  barcode.ids <- unique(TRA.data$barcode)
  
  scores <- vapply(barcode.ids, function(barcode) {
    
    TRA.subset <- TRA.data[TRA.data$barcode == barcode,]
    TRB.subset <- TRB.data[TRB.data$barcode == barcode,]
    
    if (nrow(TRA.subset) == 0 | nrow(TRB.subset) == 0) return(0)
    
    TRA.v <- TRA.subset$v
    TRA.j <- TRA.subset$j
    TRB.v <- TRB.subset$v
    TRA.length <- nchar(TRA.subset$cdr3_aa)
    
    if (is.null(TRA.v) | is.null(TRA.j) | is.null(TRB.v) | is.null(TRA.length)) return(0)
    
    TRA.v.match <- any(TRA.v %in% species.criteria$TRA.v)
    TRA.j.match <- any(TRA.j %in% species.criteria$TRA.j)
    TRA.length.match <- any(TRA.length %in% species.criteria$TRA.length)
    TRB.v.match <- ifelse(is.null(species.criteria$TRB.v), TRUE, any(TRB.v %in% species.criteria$TRB.v))
    
    as.integer(TRA.v.match & TRA.j.match & TRA.length.match & TRB.v.match)
    
  }, integer(1))
  
  output <- data.frame(barcode = barcode.ids, score = scores)
  colnames(output)[2] <- paste0(type, ".score")
  
  return(output)
}

# Internal criteria definitions
.MAIT.criteria <- list(
  mouse = list(TRA.v = "TRAV1", 
               TRA.j = "TRAJ33", 
               TRA.length = 12, 
               TRB.v = c("TRBV13", "TRBV19")), 
  human = list(TRA.v = "TRAV1-2", 
               TRA.j = c("TRAJ33", "TRAJ20", "TRAJ12"), 
               TRA.length = 12,
               TRB.v = c("TRBV6", "TRBV20"))
)

.iNKT.criteria <- list(
  mouse = list(TRA.v = "TRAV11", 
               TRA.j = "TRAJ18", 
               TRA.length = 15,
               TRB.v = NULL),
  human = list(TRA.v = "TRAV10", 
               TRA.j = "TRAJ18",
               TRA.length = c(14, 15, 16),
               TRB.v = "TRBV25-1")
)
