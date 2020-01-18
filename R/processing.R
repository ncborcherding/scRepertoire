#contigs is the individual output from cell ranger
#column refers to the default column for barcodes
#connect is the type of character in which is attaching the defualt barcode with any other characters
#' @export
stripBarcode <- function(contigs, column = 1, connector = "_", num_connects = 3) {
    count <- as.data.frame(t(data.frame(strsplit(contigs[,column], paste("['", connector, "']", sep="")), stringsAsFactors = FALSE)), stringsAsFactors = FALSE)[num_connects]
    contigs[,column] <- count
    return(contigs)
}
