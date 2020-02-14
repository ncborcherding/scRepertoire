#' Removing any additional prefixes to the barcodes of the filtered contigs.
#'
#' @param contigs The raw loaded filtered_contig_annotation.csv
#' @param column The column in which the barcodes are listed
#' @param connector The type of character in which is attaching the defualt barcode with any other characters
#' @param num_of_connects The number of strings combined with the connectors

#' @export
stripBarcode <- function(contigs, column = 1, connector = "_", num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], paste("['", connector, "']", sep="")), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#' Adding variables after the combination of contigs.
#'
#' @param df The product of CombineContig()
#' @param name The column header you'd like to add
#' @param variables The exact values to add to each element of the list

#' @export
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

#' Adding variables aftering the combination of contigs.
#'
#' @param df The product of CombineContig()
#' @param name The column header you'd like to use to subset
#' @param variables The values to subset by

#' @export
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

#This suppressing outputs for using dput()
#' @export
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

#' Allows users to take the meta data in seurat and place it into a list that will work with all the functions
#' @param df seurat object after combineSeurat()
#'
#' @export
seurat2List <- function(df) {
    if (class(df)[1] != "Seurat") {
        stop("Use a seurat object to convert into a list")
    }
    meta <- data.frame(df@meta.data, df@active.ident)
    colnames(meta)[length(meta)] <- "cluster"
    unique <- stringr::str_sort(as.character(unique(meta$cluster)), numeric = TRUE)
    df <- NULL
    for (i in seq_along(unique)) {
        subset <- subset(meta, meta[,"cluster"] == unique[i])
        df[[i]] <- subset
    }
    names(df) <- unique
    return(df)
}
