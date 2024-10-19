# testcases are NOT comprehensive

getTestTcrList <- function() {
    getdata("combineContigs", "combined")[1:3]
}

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
        makeRandomBcrBarcodesMatchTcr, n = doubletsPerSample
    )
}

makeRandomBcrBarcodesMatchTcr <- function(bcrDf, tcrDf, n) {
    sampleUniqueBarcodeDf(bcrDf, n) %>%
        dplyr::mutate(tcrBarcode = sampleUniqueBarcodeDf(tcrDf, n)$barcode) %>%
        dplyr::full_join(bcrDf, by = "barcode") %>%
        dplyr::mutate(
            barcode = ifelse(is.na(tcrBarcode), barcode, tcrBarcode)
        ) %>%
        dplyr::select(-tcrBarcode)
}

sampleUniqueBarcodeDf <- function(contigDf, n) {
    contigDf %>%
        dplyr::select(barcode) %>%
        dplyr::distinct() %>%
        dplyr::slice_sample(n = n)
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

    expect_true(nrow(getContigDoublets(tcr, bcr)) >= NUM_UNIQUE_DOUBLETS_PER_SAMPLE * length(tcr))
    # TODO
})
