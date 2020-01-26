#contigs is the individual output from cell ranger
#column refers to the default column for barcodes
#connect is the type of character in which is attaching the defualt barcode with any other characters
#' @export
stripBarcode <- function(contigs, column = 1, connector = "_", num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], paste("['", connector, "']", sep="")), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}

#df refers to the combined contig object after the combineContig()
#name is the name of the variable you'd like to add
#variables is the specific string to add to each list element. The length of the variables needs to match the length of the list elements
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

#df refers to the combined contig object after the combineContig()
#name is the name of the variable you'd like to subset based on
#variables is the values you'd like to isolate
#' @export
subsetContig <- function(df,
                        name,
                        variables = NULL) {
    df2 <- list()
    for (i in seq_along(df)) {
        if (df[[i]][1,name] %in% values) {
            df2 <- append(df2, list(df[[i]]))
        }
        else {
            next()
        }
    }
    return(df2)
}
