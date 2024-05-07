#' Generate data frame to be used with circlize R package to visualize
#' clones as a chord diagram. 
#' 
#' This function will take the meta data from the product of 
#' \code{\link{combineExpression}} and generate a relational data frame to 
#' be used for a chord diagram. Each cord will represent the number of 
#' clone unique and shared across the multiple \strong{group.by} variable. 
#' If using the downstream circlize R package, please read and cite the
#' following \href{https://pubmed.ncbi.nlm.nih.gov/24930139/}{manuscript}.
#' If looking for more advance ways for circular visualizations, there
#' is a great \href{https://jokergoo.github.io/circlize_book/book/}{cookbook}
#' for the circlize package.
#' 
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined, 
#'                                    scRep_example)
#' 
#' #Getting data frame output for Circlize
#' circles <- getCirclize(scRep_example, 
#'                        group.by = "seurat_clusters")
#' 
#' 
#' @param sc.data The single-cell object after \code{\link{combineExpression}}.
#' @param cloneCall How to call the clone - VDJC gene (\strong{gene}), 
#' CDR3 nucleotide (\strong{nt}), CDR3 amino acid (\strong{aa}),
#' VDJC gene + CDR3 nucleotide (\strong{strict}) or a custom variable 
#' in the data.  
#' @param group.by The group header for which you would like to analyze 
#' the data.
#' @param proportion Calculate the relationship unique 
#' clones (proportion = FALSE) or normalized by 
#' proportion (proportion = TRUE)
#' @param include.self Include counting the clones within a single group.by
#' comparison
#' 
#' @export
#' @concept SC_Functions
#' @return A data frame of shared clones between groups formated for \link[circlize]{chordDiagram}
#' @author Dillon Corvino, Nick Borcherding
getCirclize <- function(sc.data, 
                        cloneCall = "strict", 
                        group.by = NULL, 
                        proportion = FALSE,
                        include.self = TRUE) {
  meta <- .grabMeta(sc.data)
  cloneCall <- .theCall(meta, cloneCall)
  if(is.null(group.by)) {
    group.by <- "ident"
  }
  
  #Making exhaustive group.by dat frame
  group_pairs <- expand.grid(group1 = unique(meta[,group.by]), group2 = unique(meta[,group.by]))
  group_pairs <- unique(t(apply(group_pairs, 1, function(x) str_sort(x, numeric = TRUE))))
  group_pairs <- as.data.frame(group_pairs)
  colnames(group_pairs) <- c("from", "to")
  
  if(!include.self) {
    group_pairs <- group_pairs[group_pairs[,1] != group_pairs[,2],]
  }
  
  #Count clones across all identities
  clone.table <- .clone.counter(meta, group.by, cloneCall)
  
  group_pairs$value <- NA
  
  for(i in seq_len(nrow(group_pairs))) {
    pair1 <- group_pairs[i,1]
    pair2 <- group_pairs[i,2]
    
    clone1 <- clone.table[clone.table[,1] == pair1, cloneCall]
    clone2 <- clone.table[clone.table[,1] == pair2, cloneCall]
    
    common <- intersect(clone1, clone2)
    value <- length(common)
    
    if(pair1 == pair2) {
      tmp <- clone.table[clone.table[,cloneCall] %in% common & clone.table[,1] != pair1,]
      shared.clones <- unique(tmp[,cloneCall])
      value <- value - length(shared.clones)
    }
    
    if (proportion) {
      value <- value/length(unique(clone.table[clone.table[,1] == pair2,cloneCall]))
    }
    
    group_pairs$value[i] <- value
  }
  return(group_pairs)
}
