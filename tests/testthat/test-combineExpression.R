# test script for combineExpression.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)


test_that("Function correctly adds metadata to Seurat object", {
  sc_new <- combineExpression(combined, scRep_example, cloneCall = "strict")
  expect_s4_class(sc_new, "Seurat")
  new_meta <- sc_new@meta.data
  expected_cols <- c("CTstrict", "clonalProportion", "clonalFrequency", "cloneSize")
  expect_true(all(expected_cols %in% names(new_meta)))
  cell1_meta <- new_meta["P17B_AAACCTGCAACACGCC-1", ]
  expect_equal(cell1_meta$CTaa, "CAYRSAQAGGTSYGKLTF_CAISEQGKGELFF")
  expect_equal(cell1_meta$clonalFrequency, 1)
  expect_equal(cell1_meta$clonalProportion, 0.02, tolerance = 1e-1) 
})

test_that("`group.by` calculates frequency across specified groups", {
  # Add 'Type' to the Seurat object metadata to match the contig list
  scRep_example$Type <- c(rep("T1", 4), rep("T2", 4), "T1")
  
  sc_new <- combineExpression(combined, scRep_example, cloneCall = "strict", group.by = "Type")
  
  # Check frequency for clone 'A' which is in both T1 and T2
  # In group T1, clone A appears twice.
  cell1_meta <- sc_new@meta.data["S1_C1",]
  expect_equal(cell1_meta$clonalFrequency, 2)
  
  # In group T2, clone A appears once.
  cell5_meta <- sc_new@meta.data["S2_C5",]
  expect_equal(cell5_meta$clonalFrequency, 1)
})

test_that("`proportion = FALSE` uses frequency for `cloneSize` binning", {
  sc_new <- combineExpression(combined, scRep_example, cloneCall = "strict", proportion = FALSE,
                              cloneSize = c(Small = 1, Medium = 2, Large = 3))
  expectEqual(levels(sc_new$cloneSize), c("Large (2 < X <= 11)", 
                                          "Medium (1 < X <= 2)", "Small (0 < X <= 1)",  
                                          "None ( < X <= 0)"))
})



test_that("`filterNA = TRUE` correctly subsets the Seurat object", {
  sc_filtered <- combineExpression(combined, scRep_example, cloneCall = "strict", filterNA = TRUE)
  expect_equal(ncol(sc_filtered), 365)
})

