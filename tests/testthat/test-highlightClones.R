# test script for highlightClones.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)

test_that("Function correctly highlights a single sequence", {
  # Highlight cloneB
  sc_highlighted <- highlightClones(scRep_example, cloneCall = "aa", sequence = "CVVSDNTGGFKTIF_CASSVRRERANTGELFF")
  
  # Check that output is a Seurat object and has the new column
  expect_s4_class(sc_highlighted, "Seurat")
  expect_true("highlight" %in% names(sc_highlighted@meta.data))
  
  # Verify the contents of the 'highlight' column
  highlight_col <- sc_highlighted$highlight
  cloneB_indices <- which(scRep_example$CTaa == "CVVSDNTGGFKTIF_CASSVRRERANTGELFF")
  other_indices <- which(scRep_example$CTaa != "CVVSDNTGGFKTIF_CASSVRRERANTGELFF")
  
  expect_equal(unname(highlight_col[cloneB_indices]), rep("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", 11))
})



