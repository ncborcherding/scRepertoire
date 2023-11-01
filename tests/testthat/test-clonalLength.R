# test script for clonalLength.R - testcases are NOT comprehensive!

# testdata: (assumes combineTCR works)
combined <- combineTCR(contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))


test_that("clonalLength works", {
  expect_doppelganger(
    "clonalLength_both_chain_plot", clonalLength(combined, chain = "both") 
  )
})

