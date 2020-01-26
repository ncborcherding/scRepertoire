#' Removing any additional prefixes to the barcodes of the filtered contigs.
#'
#' @param contigs raw loaded filtered_contig_annotation.csv
#' @param column is column in which the barcodes are listed
#' @param connector is the type of character in which is attaching the defualt barcode with any other characters
#' @param num_of_connects is the number of strings combined with the connectors
#' @example stripBarcode(csv1, column = 1, connector = "_", num_connects=3)

#' @export
stripBarcode <- function(contigs, column = 1, connector = "_", num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], paste("['", connector, "']", sep="")), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#' Adding variables aftering the combination of contigs
#'
#' @param df The product of CombineContig()
#' @param name is the column header you'd like to add
#' @param variables is the exact values to add to each element of the list
#' @example addVariable(combined, name="batch", variable=c("b1","b1","b2","b2","b2","b2"))

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

#' Adding variables aftering the combination of contigs
#'
#' @param df The product of CombineContig()
#' @param name is the column header you'd like to use to subset
#' @param variables is the values to subset by
#' @example subsetContig(combined, name="sample", variable="P1")

#' @export
subsetContig <- function(df,
                        name,
                        variables = NULL) {
    df2 <- list()
    for (i in seq_along(df)) {
        if (df[[i]][1,name] %in% variables) {
            df2 <- append(df2, list(df[[i]]))
        }
        else {
            next()
        }
    }
    return(df2)
}
