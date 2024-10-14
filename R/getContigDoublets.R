#' Get Contig Doublets
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' 
#' This function identifies potential doublets by finding common barcodes
#' between TCR and BCR outputs. It extracts unique barcodes from each list
#' of dataframes, finds the intersection of the barcodes, and joins the
#' resulting data.
#'
#' @param tcrOutput Output of [combineTCR()]. A list of data.frames containing TCR contig
#' information, each dataframe must have a `barcode` column.
#' @param bcrOutput Output of [combineBCR()]. A list of data.frames containing BCR contig
#' information, each dataframe must have a `barcode` column.
#'
#' @return A dataframe of barcodes that exist in both the TCR and BCR data, with
#' columns from both sets of data. If there are no doublets, the returned
#' data.frame will have the same colnames but no rows.
#' @export
#' @examples
#' # TODO
getContigDoublets <- function(tcrOutput, bcrOutput) {

    assert_that(isListOfNonEmptyDataFrames(tcrOutput))
    assert_that(isListOfNonEmptyDataFrames(bcrOutput))

    rawBarcodeColname <- tempColnameForDfList(
        c(tcrOutput, bcrOutput), "raw_barcode"
    )

    listOfTcrBcrWithRawBarcode <- list(tcrOutput, bcrOutput) %>%
        lapplyOnAll(function(df) {
            df[[rawBarcodeColname]] <- extractBarcodeStrings(df$barcode)
            df
        })

    doubletBarcodes <- listOfTcrBcrWithRawBarcode %>%
        lapplyOnAll(function(df) {
            unique(df[[rawBarcodeColname]])
        }) %>%
        lapply(purrr::list_flatten) %>%
        purrr::reduce(intersect)

    if (length(doubletBarcodes) == 0) {
        return(makeEmptyIntersectionDf(tcrOutput[[1]], bcrOutput[[1]]))
    }

    listOfTcrBcrWithRawBarcode %>%
        lapplyOnAll(function(df) {
            df[df[[rawBarcodeColname]] %in% doubletBarcodes, ]
        }) %>%
        lapply(dplyr::bind_rows) %>%
        purrr::reduce(autoFullJoin)
}

tempColnameForDfList <- function(dfList, baseName = "temp") {
    colnameSet <- unique(unlist(lapply(dfList, colnames)))
    tail(make.unique(c(colnameSet, baseName)), 1)
}

lapplyOnAll <- function(listOfLists, fun) {
    lapply(listOfLists, function(x) lapply(x, fun))
}

makeEmptyIntersectionDf <- function(...) {
    purrr::reduce(list(...), function(df1, df2) {
        autoFullJoin(df1[0, ], df2[0, ])
    })
}

autoFullJoin <- function(df1, df2) {
    suppressMessages(dplyr::full_join(df1, df2))
}

extractBarcodeStrings <- function(inputStrings) {
    matches <- unlist(lapply(inputStrings, function(x) {
        regmatches(x, gregexpr("[a-zA-Z_]+_[ATGC]+-\\d+", x))
    }))
    matches[matches != ""]
}
