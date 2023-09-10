#' Allows users to take the meta data in Seurat/SCE and place it into a list 
#' that will work with all the functions
#'
#' Allows users to perform more fundamental measures of clonotype analysis 
#' using the meta data from the Seurat or SCE object. For Seurat objects the 
#' active identity is automatically added as "cluster". Remaining grouping 
#' parameters or SCE or Seurat objects must appear in the meta data.
#'
#' @examples
#' #Getting the combined contigs
#' combined <- combineTCR(contig_list, 
#'                         samples = c("P17B", "P17L", "P18B", "P18L", 
#'                                     "P19B","P19L", "P20B", "P20L"))
#' 
#' #Getting a sample of a Seurat object
#' scRep_example <- get(data("scRep_example"))
#' scRep_example <- combineExpression(combined, scRep_example)
#' 
#' #Using expression2List
#' newList <- expression2List(scRep_example, split.by = "seurat_clusters")
#' 
#' @param sc object after combineExpression().
#' @param split.by The column header to group the new list. NULL will return clusters.
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of single-cell object with 
#' elements divided by the group parameter
expression2List <- function(sc, split.by) {
  if (!inherits(x=sc, what ="Seurat") & 
      !inherits(x=sc, what ="SummarizedExperiment")) {
    stop("Use a Seurat or SCE object to convert into a list")
  }
  meta <- grabMeta(sc)
  if(is.null(split.by)){
    split.by <- "cluster"
  }
  unique <- str_sort(as.character(unique(meta[,split.by])), numeric = TRUE)
  df <- NULL
  for (i in seq_along(unique)) {
    subset <- subset(meta, meta[,split.by] == unique[i])
    subset <- subset(subset, !is.na(cloneType))
    df[[i]] <- subset
  }
  names(df) <- unique
  return(df)
}