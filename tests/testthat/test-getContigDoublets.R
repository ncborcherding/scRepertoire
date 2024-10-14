test_that("getContigDoublets works for no doublets", {

    tcr <- getdata("combineContigs", "combineTCR_list_expected")

    # create a BCR list from testdata with no doublets
    bcr <- getdata("combineContigs", "combineBCR_list_expected")
    bcr <- list(bcr[[1]][1:10, ], bcr[[1]][20:30, ], bcr[[1]][100:110, ])
    names(bcr) <- names(tcr)

    expected_no_doublet_output <- structure(
        list(
            barcode = character(0), sample = character(0),
            TCR1 = character(0), cdr3_aa1 = character(0),
            cdr3_nt1 = character(0), TCR2 = character(0),
            cdr3_aa2 = character(0), cdr3_nt2 = character(0),
            CTgene = character(0), CTnt = character(0), CTaa = character(0),
            CTstrict = character(0), IGH = character(0), IGLC = character(0),
            contigType = structure(
                integer(0), levels = c("BCR", "TCR"), class = "factor"
            )
        ),
        row.names = integer(0),
        class = "data.frame"
    )

    expect_identical(getContigDoublets(tcr, bcr), expected_no_doublet_output)
})

test_that("getContigDoublets works for inputs with doublets", {

    tcr <- getdata("combineContigs", "combineTCR_list_expected")

    # create a BCR list from testdata with doublets
    bcr <- getdata("combineContigs", "combineBCR_list_expected")
    bcr <- list(bcr[[1]][1:10, ], bcr[[1]][20:30, ], bcr[[1]][100:110, ])
    names(bcr) <- names(tcr)

    # UNFINISHED
    # TODO purposely introduce doublets into testing data
    # TODO test expected output
    expect_equal(2 * 2, 4)
})
