# test script for combineExpression.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)


test_that("Basic functionality with Seurat object and default parameters", {
  # Run the function with default settings
  sc <- combineExpression(combined, scRep_example)
  expect_s4_class(sc, "Seurat")
  expect_equal(ncol(sc), ncol(scRep_example))
  metadata_cols <- colnames(sc[[]])
  expected_cols <- c("CTstrict", "clonalFrequency", "clonalProportion", "cloneSize")
  expect_true(all(expected_cols %in% metadata_cols))
  expect_s3_class(sc$cloneSize, "factor")
})


test_that("`cloneCall` parameter works correctly", {
  # Test with cloneCall = "nt"
  sc_nt <- combineExpression(combined, scRep_example, cloneCall = "nt")
  expect_true("CTnt" %in% colnames(sc_nt[[]]))
  
  # Test with cloneCall = "aa"
  sc_aa <- combineExpression(combined, scRep_example, cloneCall = "aa")
  expect_true("CTaa" %in% colnames(sc_aa[[]]))
  
  # Test with cloneCall = "gene"
  sc_gene <- combineExpression(combined, scRep_example, cloneCall = "gene")
  expect_true("CTgene" %in% colnames(sc_gene[[]]))
})


test_that("`proportion` and `cloneSize` parameters work as expected", {
  # Test with frequency calculation (proportion = FALSE)
  sc_freq <- combineExpression(combined, 
                               scRep_example, 
                               proportion = FALSE, 
                               cloneSize = c(1,2,5,10))
  meta_freq <- sc_freq[[]]
  max_freq <- max(meta_freq$clonalFrequency, na.rm = TRUE)
  clone_id_max_freq <- meta_freq$CTstrict[which.max(meta_freq$clonalFrequency)]
  
  # Default cloneSize for frequency goes up to the max value automatically
  expect_true(grepl(paste0("<= ", max_freq, ")"), levels(meta_freq$cloneSize)[1]))
  
  # Test error when proportion = FALSE and cloneSize bins are < 1
  bad_cloneSize <- c(Rare = 0.1, Small = 0.5)
  expect_error(
    combineExpression(combined, scRep_example, proportion = FALSE, cloneSize = bad_cloneSize),
    "Adjust the cloneSize parameter - there are groupings < 1"
  )
})

test_that("`filterNA` parameter correctly subsets the object", {
  sc_filtered <- combineExpression(combined, scRep_example, filterNA = TRUE)
  expect_lt(ncol(sc_filtered), ncol(scRep_example))
  expect_false(any(is.na(sc_filtered$CTstrict)))
  non_na_barcodes <- rownames(scRep_example[[]])[!is.na(combineExpression(combined, scRep_example)$CTstrict)]
  expect_equal(ncol(sc_filtered), length(non_na_barcodes))
})

test_that("`chain` parameter filters correctly for calculations", {
  sc_tra <- combineExpression(combined, scRep_example, chain = "TRA")
  sc_trb <- combineExpression(combined, scRep_example, chain = "TRB")
  
  # The clonal proportions should be different when calculated on different chains
  expect_false(identical(sc_tra$clonalProportion, sc_trb$clonalProportion))
  
  # The final clone ID (CTstrict) should still represent the full paired chain
  # This is a key feature of the function
  sc_both <- combineExpression(combined, scRep_example, chain = "both")
  
  # Pick a cell to compare
  cell_barcode <- na.omit(sc_tra[[]])$barcode[1]
  
  # But the strict clonotype ID should be identical, preserving the paired info
  id_tra <- sc_tra[[c("CTstrict")]][cell_barcode, ]
  id_both <- sc_both[[c("CTstrict")]][cell_barcode, ]
  expect_equal(id_tra, id_both)
})

test_that("Function issues a warning for high barcode mismatch", {
  # Create a Seurat object with non-matching barcodes
  combined.mm <- lapply(combined, function(x) {
                x$barcode <- paste0("mismatch_", x$barcode )
                x
  })
  
  # Expect the specific warning message
  expect_warning(
    combineExpression(combined.mm, scRep_example),
    "< 1% of barcodes match"
  )
})
