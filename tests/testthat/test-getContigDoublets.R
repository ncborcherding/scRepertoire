# testcases are NOT comprehensive

getTestTcrList <- function() {
    getdata("combineContigs", "combined")[1:3]
}

#' Generate a [combineBCR()] output example with matching sample names
#' as the output of [getTestTcrList()] above and no doublets.
getTestBcrListNoDoublets <- function() {
    bcr <- getdata("combineContigs", "combineBCR_list_expected")
    bcr <- list(bcr[[1]][1:10, ], bcr[[1]][20:30, ], bcr[[1]][100:110, ])
    names(bcr) <- names(getTestTcrList())
    purrr::imap(bcr, function(df, sampleName) {
        df$sample <- sampleName
        df
    })
}

getTestBcrListWithDoublets <- function(doubletsPerSample, seed = 42) {
    if (!is.null(seed)) withr::local_seed(seed)
    purrr:::map2(
        getTestBcrListNoDoublets(), getTestTcrList(),
        makeRandomBcrBarcodesMatchTcr,
        n = doubletsPerSample
    )
}

makeRandomBcrBarcodesMatchTcr <- function(bcrDf, tcrDf, n) {
    sampleUniqueBarcodeDf(bcrDf, n) %>%
        dplyr::mutate(
            tcrBarcode = sampleUniqueBarcodeDf(tcrDf, n, asDf = FALSE)
        ) %>%
        dplyr::full_join(bcrDf, by = "barcode") %>%
        dplyr::mutate(
            barcode = ifelse(is.na(tcrBarcode), barcode, tcrBarcode)
        ) %>%
        dplyr::select(-tcrBarcode)
}

sampleUniqueBarcodeDf <- function(contigDf, n, asDf = TRUE) {
    contigDf %>%
        dplyr::select(barcode) %>%
        dplyr::distinct() %>%
        dplyr::slice_sample(n = n) %>%
        (if (asDf) identity else function(x) x$barcode)
}

test_that("getContigDoublets works for no doublets", {

    expected_no_doublet_output <- structure(
        list(
            contigType = structure(
                integer(0), levels = c("BCR", "TCR"), class = "factor"
            ),
            barcode = character(0), sample = character(0),
            TCR1 = character(0), cdr3_aa1 = character(0),
            cdr3_nt1 = character(0), TCR2 = character(0),
            cdr3_aa2 = character(0), cdr3_nt2 = character(0),
            CTgene = character(0), CTnt = character(0), CTaa = character(0),
            CTstrict = character(0), IGH = character(0), IGLC = character(0)
        ),
        row.names = integer(0),
        class = "data.frame"
    )

    expect_identical(
        getContigDoublets(getTestTcrList(), getTestBcrListNoDoublets()),
        expected_no_doublet_output
    )
})

test_that("getContigDoublets works for inputs with doublets", {

    NUM_UNIQUE_DOUBLETS_PER_SAMPLE <- 3
    tcr <- getTestTcrList()
    bcr <- getTestBcrListWithDoublets(NUM_UNIQUE_DOUBLETS_PER_SAMPLE)

    doubletDf <- getContigDoublets(tcr, bcr)

    expect_equal(
        nrow(doubletDf),
        NUM_UNIQUE_DOUBLETS_PER_SAMPLE * length(tcr) * 2
    )

    expect_identical(
        colnames(doubletDf),
        c("contigType", "barcode", "sample", "TCR1", "cdr3_aa1", "cdr3_nt1",
          "TCR2", "cdr3_aa2", "cdr3_nt2", "CTgene", "CTnt", "CTaa", "CTstrict",
          "IGH", "IGLC")
    )

    getBarcodeSampleForContigType <- function(contigType) {
        doubletDf %>%
            dplyr::filter(contigType == contigType) %>%
            dplyr::select(barcode, sample) %>%
            dplyr::arrange(barcode, sample)
    }

    expect_identical(
        getBarcodeSampleForContigType("BCR"),
        getBarcodeSampleForContigType("TCR")
    )

    makeCharNaDf <- function(dfColnames, nrow) {
        matrix(nrow = nrow, ncol = length(dfColnames)) %>%
            data.frame() %>%
            (function(df) {
                colnames(df) <- dfColnames
                df
            }) %>%
            dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    }

    expect_identical(
        doubletDf %>%
            dplyr::filter(contigType == "TCR") %>%
            dplyr::select(IGH, IGLC),
        makeCharNaDf(
            c("IGH", "IGLC"), NUM_UNIQUE_DOUBLETS_PER_SAMPLE * length(tcr)
        )
    )

    expect_identical(
        doubletDf %>%
            dplyr::filter(contigType == "BCR") %>%
            dplyr::select(TCR1, TCR2),
        makeCharNaDf(
            c("TCR1", "TCR2"),
            NUM_UNIQUE_DOUBLETS_PER_SAMPLE * length(bcr)
        )
    )

})
