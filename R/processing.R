#' Removing any additional prefixes to the barcodes of the filtered contigs.
#'
#' @param contigs The raw loaded filtered_contig_annotation.csv
#' @param column The column in which the barcodes are listed
#' @param connector The type of character in which is attaching the defualt barcode with any other characters
#' @param num_connects The number of strings combined with the connectors

#' @export
#' @return list with the suffixes of the barcodes removed.
stripBarcode <- function(contigs, column = 1, connector = "_", num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], paste("['", connector, "']", sep="")), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#' Adding variables after the combination of contigs.
#'
#' This function adds variables to the product of combineContig() to be used in later visualizations. For each element
#' in the combinContig(), the function will add a column (labled by name) with the variable. The length of the
#' variable paramater needs to match the length of the combineContig() object.
#'
#' @examples
#' combined <- addVariable(combined, name = "batch", variables = c(1,1,2,3,3))
#'
#' @param df The product of CombineContig().
#' @param name The column header to add.
#' @param variables The exact values to add to each element of the list.

#' @export
#' @return list of contigs with a new column (name).
addVariable <- function(df,
                        name = NULL,
                        variables =  NULL) {
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
#' This function allows for the subsetting of the product of combineContig() by the name of the individual list element.
#' In gneral the names of are samples + _ + ID, allowing for users to subset the product of combineContig() across a
#' string or indivudal name.
#'
#' @examples
#' subset <- subsetContig(combined, name = "sample", variables = c("X", "Y")
#'
#' @param df The product of CombineContig().
#' @param name The column header you'd like to use to subset.
#' @param variables The values to subset by, must be in the names(df).

#' @export
#' @return list of contigs that have been filtered for the name parameter
subsetContig <- function(df,
                        name,
                        variables = NULL) {
    names2 <- NULL
    df2 <- list()
    for (i in seq_along(df)) {
        if (df[[i]][1,name] %in% variables) {
            df2 <- append(df2, list(df[[i]]))
            n2 <- names(df)[i] #need to store the name of the list element and replace
            names2 <- c(names2, n2)
        }
        else {
            next()
        }
    }
    names(df2) <- names2
    return(df2)
}

# This suppressing outputs for using dput()
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

#This is to help sort the type of clonotype data to use
theCall <- function(x) {
if (x == "gene") {
    x <- "CTgene"
} else if(x == "nt") {
    x <- "CTnt"
} else if (x == "aa") {
    x <- "CTaa"
} else if (x == "gene+nt") {
    x <- "CTstrict"
}
    return(x)
}

#' Allows users to take the meta data in seurat/SCE and place it into a list that will work with all the functions
#'
#' Allows users to perform more funadmental measures of clonotype analysis using the meta data from the seurat or SCE
#' object. For Seurat objects the active identity is automatically added as "cluster". Reamining grouping parameters
#' for SCE or Seurat objects must appear in the meta data.
#'
#' @examples
#' newList <- expression2List(seurat, group = "cluster")
#'
#' @param sc object after combineExpression().
#' @param group The column header to group the new list by
#' @importFrom stringr str_sort
#' @export
#' @return list derived from the meta data of single-cell object with elements divided by the group parameter
expression2List <- function(sc, group) {
    if (!inherits(x=sc, what ="Seurat") & !inherits(x=sc, what ="SummarizedExperiment")) {
        stop("Use a seurat or SCE object to convert into a list")
    }

    if (inherits(x=sc, what ="Seurat")) {
        meta <- data.frame(sc[[]], Idents(sc))
        colnames(meta)[length(meta)] <- "cluster"
    }
    else if (inherits(x=sc, what ="SummarizedExperiment")){
        meta <- sc@metadata[[1]]
    }
    unique <- str_sort(as.character(unique(meta[,group])), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,"cluster"] == unique[i])
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}
