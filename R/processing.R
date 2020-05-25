#' Removing any additional prefixes to the barcodes of filtered contigs.
#'
#' @param contigs The raw loaded filtered_contig_annotation.csv
#' @param column The column in which the barcodes are listed
#' @param connector The type of character in which is attaching the defualt 
#' barcode with any other characters
#' @param num_connects The number of strings combined with the connectors
#' @examples 
#' stripBarcode(contig_list[[1]], column = 1, connector = "_", num_connects = 3)
#' @export
#' @return list with the suffixes of the barcodes removed.
stripBarcode <- function(contigs, column = 1, connector = "_", 
                            num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], 
                            paste("['", connector, "']", sep="")), 
                            stringsAsFactors = FALSE)), 
                            stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#' Adding variables after the combination of contigs.
#'
#' This function adds variables to the product of combineContig() to be used in 
#' later visualizations. For each element in the combinContig(), the function 
#' will add a column (labled by name) with the variable. The length of the
#' variable paramater needs to match the length of the combineContig() object.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' combined <- addVariable(combined, name = "batch", variables = c(1,1,1,1,2,2))
#'
#' @param df The product of CombineContig().
#' @param name The column header to add.
#' @param variables The exact values to add to each element of the list.
#' @export
#' @return list of contigs with a new column (name).
addVariable <- function(df, name = NULL, variables =  NULL) {
    if (length(df) != length(variables)) {
        stop("Make sure the variables match the length of the contig list")
    }
    for (i in seq_along(df)) {
        df[[i]][,name] <- variables[i]
    }
    return(df)
}

#' Subset the combineContig() product.
#'
#' This function allows for the subsetting of the product of combineContig() by 
#' the name of the individual list element. In general the names of are samples 
#' + _ + ID, allowing for users to subset the product of combineContig() across 
#' a string or individual name.
#'
#' @examples
#' x <- contig_list
#' combined <- combineTCR(x, rep(c("PX", "PY", "PZ"), each=2), 
#' rep(c("P", "T"), 3), cells ="T-AB")
#' subset <- subsetContig(combined, name = "sample", variables = c("PX"))
#'
#' @param df The product of CombineContig().
#' @param name The column header you'd like to use to subset.
#' @param variables The values to subset by, must be in the names(df).

#' @export
#' @return list of contigs that have been filtered for the name parameter
subsetContig <- function(df, name, variables = NULL) {
    names2 <- NULL
    df2 <- list()
    for (i in seq_along(df)) {
        if (df[[i]][1,name] %in% variables) {
            df2 <- append(df2, list(df[[i]]))
            n2 <- names(df)[i] 
            names2 <- c(names2, n2)
        }
        else {
            next()
        }
    }
    names(df2) <- names2
    return(df2)
}

#' Allows users to take the meta data in seurat/SCE and place it into a list 
#' that will work with all the functions
#'
#' Allows users to perform more fundamental measures of clonotype analysis 
#' using the meta data from the seurat or SCE object. For Seurat objects the 
#' active identity is automatically added as "cluster". Reamining grouping 
#' parameters or SCE or Seurat objects must appear in the meta data.
#'
#' @examples
#' \donttest{
#' newList <- expression2List(seurat, group = "cluster")
#' }
#' @param sc object after combineExpression().
#' @param group The column header to group the new list by
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of single-cell object with 
#' elements divided by the group parameter
expression2List <- function(sc, group) {
    if (!inherits(x=sc, what ="Seurat") & 
        !inherits(x=sc, what ="SummarizedExperiment")) {
            stop("Use a seurat or SCE object to convert into a list")
    }
    meta <- grabMeta(sc)
    unique <- str_sort(as.character(unique(meta[,group])), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,"cluster"] == unique[i])
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}

