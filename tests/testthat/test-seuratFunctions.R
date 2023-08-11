# test script for seuratFunctions.R - testcases are NOT comprehensive!

test_that("combineExpression works with seurat objects", {
	data("combined_mini_contig_list", "screp_example")
	combined_test <- combineExpression(combined_mini_contig_list, screp_example)
	
	expect_length(combined_test@meta.data, 15)
	expect_equal(combined_test@meta.data[, 1:8], screp_example@meta.data[, 1:8])
	expect_equal(
		combined_test@meta.data[, 9:15],
		getdata("seuratFunctions", "combineExpression_new_metadata")
	)
})

# TODO more testcases for combineEXpression, especially with SCE objects
# TODO highlightClonotypes
# TODO alluvialClonotypes
# TODO occupiedscRepertoire
# TODO clonalOverlay
# TODO createHTOContigList
