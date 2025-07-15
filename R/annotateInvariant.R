#' Annotate invariant T cells (MAIT or iNKT) in single-cell TCR data
#' 
#' The [annotateInvariant()] function identifies potential mucosal-associated 
#' invariant T (`MAIT``) cells or invariant natural killer T (`iNKT`) cells from 
#' single-cell sequencing datasets based on their characteristic  TCR usage. 
#' It extracts TCR chain information from the provided single-cell 
#' data, checks it against known invariant T-cell receptor criteria for either 
#' MAIT or iNKT cells, and returns a score indicating the presence (1) or 
#' absence (0) of these invariant cell populations for each individual cell. 
#' The function supports data from mouse and human samples, providing a 
#' convenient method to annotate specialized T-cell subsets within single-cell 
#' analyses.
#'
#' @param input.data The product of [combineTCR()] or [combineExpression()].
#' @param type Character specifying the type of invariant cells to
#'annotate (`MAIT` or `iNKT`).
#' @param species Character specifying the species 
#' ('mouse' or 'human').
#' 
#' @examples
#' # Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' # Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' 
#' # Using combineExpresion()
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' # Using annotateInvariant()
#' annotateInvariant(input.data = scRep_example, type = "MAIT", species = "human")
#' annotateInvariant(input.data = scRep_example, type = "iNKT", species = "human")
#' 
#' @importFrom immApex getIR
#' @importFrom  rlang %||%
#' @export
#' @return A single-cell object or list with the corresponding annotation 
#' scores (0 or 1) added.
annotateInvariant <- function(input.data, 
                              type = c("MAIT", "iNKT"), 
                              species = c("mouse", "human")) {
  
  type <- match.arg(type)
  species <- match.arg(species)
  
  if(!.is.seurat.or.se.object(input.data) & !inherits(input.data, "list")) {
    stop("Please use the output of combineTCR() or combineExpression() as input.data")
  }

  
  TCRS <- lapply(c("TRA", "TRB"), function(x) {
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
    
    TRA.v.match <- grepl(paste(species.criteria$TRA.v, collapse = "|"), TRA.v)
    TRA.j.match <- grepl(paste(species.criteria$TRA.j, collapse = "|"), TRA.j)
    TRA.length.match <- any(TRA.length %in% species.criteria$TRA.length)
    TRB.v.match <- ifelse(is.null(species.criteria$TRB.v), TRUE, grepl(paste(species.criteria$TRB.v, collapse = "|"), TRB.v))
    
    as.integer(TRA.v.match & TRA.j.match & TRA.length.match & TRB.v.match)
    
  }, integer(1))
  
  output <- data.frame(row.names = barcode.ids, score = scores)
  new.variable.name <- paste0(type, ".score")
  colnames(output)[1] <- new.variable.name
  
  if(.is.seurat.or.se.object(input.data)) {
    if (.is.seurat.object(input.data)) { 
      input.data[[new.variable.name]] <- output
    } else {
      combined_col_names <- unique(c(colnames(colData(input.data)), new.variable.name))
      full_data <- merge(colData(input.data), output, by = "row.names", all.x = TRUE)
      # at this point, the rows in full_data are shuffled. match back with the original colData
      full_data <- full_data[match(colnames(input.data), full_data[,1]), ]
      rownames(full_data) <- full_data[, 1]
      full_data  <- full_data[, -1]
      colData(input.data) <- DataFrame(full_data[, combined_col_names])  
    }
  } else {
    input.data <- lapply(input.data, function(x) {
      merged_x <- merge(x, output, by.x = "barcode", by.y = "row.names", all.x = TRUE)  
      return(merged_x)
    })
  }
  
  return(input.data)
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
