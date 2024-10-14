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
#' @param tcrOutput Output of [combineTCR()]. A list of data.frames containing
#' TCR contig information, each dataframe must have a `barcode` column.
#' @param bcrOutput Output of [combineBCR()]. A list of data.frames containing
#' BCR contig information, each dataframe must have a `barcode` column.
#'
#' @return
#' A dataframe of barcodes that exist in both the TCR and BCR data, with
#' columns from both sets of data. There will be an additional column
#' `contigType` of type factor with levels 'TCR' and 'BCR' indicating the
#' origin of the contig - this will be the new first column.
#'
#' If there are no doublets, the returned
#' data.frame will have the same colnames but no rows.
#'
#' @export
getContigDoublets <- function(tcrOutput, bcrOutput) {

    assert_that(isListOfNonEmptyDataFrames(tcrOutput))
    assert_that(all(sapply(tcrOutput, function(d) "barcode" %in% colnames(d))))
    assert_that(isListOfNonEmptyDataFrames(bcrOutput))
    assert_that(all(sapply(bcrOutput, function(d) "barcode" %in% colnames(d))))

    doubletBarcodes <- getContigDoubletBarcodes(tcrOutput, bcrOutput)

    if (length(doubletBarcodes) == 0) {
        autoFullJoin(tcrOutput[[1]][0, ], bcrOutput[[1]][0, ]) %>%
            dplyr::mutate(
                contigType = factor(character(0), levels = c("BCR", "TCR")),
                .before = 1
            ) %>%
            return()
    }

    list(TCR = tcrOutput, BCR = bcrOutput) %>%
        lapplyOnAll(function(df) {
            df[df$barcode %in% doubletBarcodes, ]
        }) %>%
        purrr::imap(function(x, type) {
            dplyr::bind_rows(x) %>%
                dplyr::mutate(contigType = type)
        }) %>%
        purrr::reduce(autoFullJoin) %>%
        dplyr::mutate(
            contigType = factor(contigType, levels = c("BCR", "TCR"))
        ) %>%
        dplyr::relocate(contigType)
}

getContigDoubletBarcodes <- function(tcrOutput, bcrOutput) {
    intersect(
        dplyr::bind_rows(tcrOutput)$barcode,
        dplyr::bind_rows(bcrOutput)$barcode
    )
}

lapplyOnAll <- function(listOfLists, fun) {
    lapply(listOfLists, function(x) lapply(x, fun))
}

autoFullJoin <- function(df1, df2) {
    suppressMessages(dplyr::full_join(df1, df2))
}
