# test script for loadContigs.R - testcases are NOT comprehensive!

test_that("loadContigs works", {
    TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
    trial1 <- loadContigs(TRUST4, format = "TRUST4")
    expect_identical(trial1, 
                     getdata("load", "loadContigs_TRUST4")
    )
    
    
    BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
    trial2 <- loadContigs(BD, format = "BD")
    expect_identical(trial2, 
                     getdata("load", "loadContigs_BD")
    )
    
    WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
    trial3 <- loadContigs(WAT3R, format = "WAT3R")
    expect_identical(trial3, 
                     getdata("load", "loadContigs_WAT3R")
    )
    
    data("contig_list")
    trial4 <- loadContigs(contig_list[[1]], format = "10X")
    expect_identical(trial4, 
                     getdata("load", "loadContigs_10x")
    )
    
    
    MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
    trial5 <- loadContigs(MIXCR, format = "MiXCR")
    expect_identical(trial5, 
                     getdata("load", "loadContigs_MiXCR")
    )
    
    Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
    trial6 <- loadContigs(Immcantation, format = "Immcantation")
    expect_identical(trial6, 
                     getdata("load", "loadContigs_Immcantation")
    )
    
    OS <- read.csv("https://www.borch.dev/uploads/contigs/OS_contigs.csv")
    trial7 <- loadContigs(OS, format = "Omniscope")
    expect_identical(trial7, 
                     getdata("load", "loadContigs_Omniscope")
    )
    
    Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
    trial8 <- loadContigs(Parse, format = "ParseBio")
    expect_identical(trial8, 
                     getdata("load", "loadContigs_Parse")
    )
}) 

# TODO Add tests for .json and AIRR and OS
# TODO Would be nice to have a dir option
