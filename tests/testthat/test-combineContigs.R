# test script for combineContigs.R - testcases are NOT comprehensive!

test_that("combineTCR works with default parameters", {
  combined <- combineTCR(contig_list[1:2], samples = c("P17B", "P17L"))
  expect_type(combined, "list")
  expect_length(combined, 2)
  expect_s3_class(combined[[1]], "data.frame")
  # Check if barcodes are prefixed
  expect_true(startsWith(combined[[1]]$barcode[1], "P17B_"))
})

test_that("combineTCR `samples` and `ID` parameters work", {
  combined <- combineTCR(contig_list[1], samples = "S1", ID = "A")
  expect_equal(names(combined)[1], "S1_A")
  expect_true(startsWith(combined[[1]]$barcode[1], "S1_A_"))
})

test_that("combineTCR `filterNonproductive = FALSE` keeps non-productive chains", {
  contig_mock <- contig_list[[1]]
  contig_mock$productive[1:50] <- "False"
  combined_filtered <- combineTCR(list(contig_mock), samples="S1")
  combined_unfiltered <- combineTCR(list(contig_mock), samples="S1", filterNonproductive = FALSE)
  expect_lt(nrow(combined_filtered[[1]]), nrow(combined_unfiltered[[1]]))
})

test_that("combineTCR `removeNA` and `removeMulti` work", {
  combined_removeNA <- combineTCR(list(contig_mock), samples="S1", removeNA = TRUE)[[1]]
  expect_true(all(!grepl("NA_", combined_removeNA$CTaa)))
  expect_true(all(!grepl("_NA", combined_removeNA$CTnt)))
  
  combined_removeMulti <- combineTCR(list(contig_mock), samples="S1", removeMulti = TRUE)
  expect_true(all(!grepl(";", combined_removeMulti$CTaa)))
  expect_true(all(!grepl(";", combined_removeMulti$CTnt)))
})


test_that("combineBCR works", {

  BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
  combined_bcr <- combineBCR(BCR, 
                    samples = "Patient1")
  expect_true(any(grepl("Cluster", combined_bcr[[1]]$CTstrict)))
  expect_type(combined_bcr, "list")
  expect_length(combined_bcr, 1)
  expect_s3_class(combined_bcr[[1]], "data.frame")
  # Check if barcodes are prefixed
  expect_true(startsWith(combined_bcr[[1]]$barcode[1], "Patient1_"))

  expect_error(combineBCR(BCR), 
               regexp = "requires the samples")

})
