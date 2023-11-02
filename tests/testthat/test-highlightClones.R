# test script for highlightClones.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
test_obj <- combineExpression(combined, scRep_example)


test_that("highlightClones works", {
  test_obj <- highlightClones(test_obj, 
                              cloneCall = "aa",
                              sequence = "CVVSDNTGGFKTIF_CASSVRRERANTGELFF")
  expect_equal(
    test_obj@meta.data[1:10, 7:14],
    getdata("seuratFunctions", "highlightClones_meta")
  )
})
