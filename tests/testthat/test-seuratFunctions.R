# test script for seuratFunctions.R - testcases are NOT comprehensive!

test_that("combineExpression works with seurat objects", {
	data("mini_contig_list", "scRep_example")
  combined <- combineTCR(mini_contig_list, 
                         samples = c("P17B", "P17L", "P18B", "P18L", 
                                     "P19B","P19L", "P20B", "P20L"))
	combined_test <- combineExpression(combined, scRep_example)
	
	expect_length(combined_test@meta.data, 12)
	expect_equal(combined_test@meta.data[, 1:6], scRep_example@meta.data[, 1:6])
	expect_equal(
		combined_test@meta.data[, 7:12],
		getdata("seuratFunctions", "combineExpression_new_metadata")
	)
})

# TODO more testcases for combineEXpression, especially with SCE objects
# TODO highlightClonotypes
# TODO alluvialClonotypes
# TODO occupiedscRepertoire
# TODO clonalOverlay
# TODO createHTOContigList
