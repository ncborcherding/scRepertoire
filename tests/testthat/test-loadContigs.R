# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works", {
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
	
}) 

# TODO Test more formats
# TODO Would be nice to have a dir option
