# test script for eportClones.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))


test_that("exportClones works", {
  trial1 <- exportClones(combined, 
                         write.file = FALSE,
                         format = "paired")
  expect_equal(
    trial1,
    getdata("export", "exportClones_paired")
  )
  
  trial2 <- exportClones(combined, 
                         write.file = FALSE,
                         format = "airr")
  expect_equal(
    trial2,
    getdata("export", "exportClones_airr")
  )
  
  trial3 <- exportClones(combined, 
                         write.file = FALSE,
                         format = "TCRMatch")
  expect_equal(
    trial3,
    getdata("export", "exportClones_TCRMatch")
  )
  
})