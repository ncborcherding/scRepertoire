# test script for eportClones.R - testcases are NOT comprehensive!

test_that("exportClones works", {

  combined <- getCombined()

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