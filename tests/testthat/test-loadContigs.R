# test script for loadContigs.R - testcases are NOT comprehensive!

test_that("loadContigs works", {
    
    BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
    trial2 <- loadContigs(BD, format = "BD")
    expect_identical(trial2, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_BD"))
    )
    
    WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
    trial3 <- loadContigs(WAT3R, format = "WAT3R")
    expect_identical(trial3, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_WAT3R"))
    )
    
    data("contig_list")
    trial4 <- loadContigs(contig_list[[1]], format = "10X")
    expect_identical(trial4, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_10x"))
    )
    
    
    MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
    trial5 <- loadContigs(MIXCR, format = "MiXCR")
    expect_identical(trial5, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_MiXCR"))
    )
    
    Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
    trial6 <- loadContigs(Immcantation, format = "Immcantation")
    expect_identical(trial6, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Immcantation"))
    )
    
    OS <- read.csv("https://www.borch.dev/uploads/contigs/OS_contigs2.csv")
    trial7 <- loadContigs(OS, format = "Omniscope")
    expect_identical(trial7, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Omniscope"))
    )
    
    Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
    trial8 <- loadContigs(Parse, format = "ParseBio")
    expect_identical(trial8, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Parse"))
    )
    
    Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
    trial9 <- loadContigs(Dandelion, format = "Dandelion")
    expect_identical(trial9, 
                     rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Dandelion"))
    )
})

test_that("loadContigs(format='TRUST4') works", {

    TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
    expect_identical(
        loadContigs(TRUST4, format = "TRUST4"), 
        rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_TRUST4"))
    )

    oneRowTrust4Input <- structure(
        list(
            `#barcode` = "CGTAGCGGTGATAAGT-1",
            cell_type = "B",
            chain1 = "*",
            chain2 = "IGKV1D-43,*,IGKJ1,IGKC,TGTCAACAGTATAGTAGTGTCCCCTGGACGTTC,CQQYSSVPWTF,6.00,CGTAGCGGTGATAAGT-1_2,76.00,0",
            secondary_chain1 = "*",
            secondary_chain2 = "*"
        ),
        row.names = c(NA, -1L),
        class = "data.frame"
    )

    expectedParsedTrust4Data <- list(
        structure(
            list(
                barcode = "CGTAGCGGTGATAAGT-1",
                v_gene = "IGKV1D-43",
                d_gene = "None",
                j_gene = "IGKJ1",
                c_gene = "IGKC",
                cdr3_nt = "TGTCAACAGTATAGTAGTGTCCCCTGGACGTTC",
                cdr3 = "CQQYSSVPWTF",
                reads = "6.00",
                chain = "IGK"
            ),
            row.names = 1L,
            class = "data.frame"
        )
    )

    expect_identical(
        loadContigs(oneRowTrust4Input, format = "TRUST4"),
        expectedParsedTrust4Data
    )
})

# TODO Add tests for .json and AIRR
# TODO Would be nice to have a dir option
