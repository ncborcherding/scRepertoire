test_that("make_screp_seurat_cmd works", {
	data("scRep_example", "mini_contig_list")
	
	test_obj <- combineExpression(mini_contig_list, scRep_example)
	expect_true(!is.null(test_obj@commands[["combineExpression"]]))
	
	test_obj <- test_obj@commands[["combineExpression"]]
	
	expect_identical(test_obj@name, "combineExpression")
	expect_identical(class(test_obj@time.stamp), c("POSIXct", "POSIXt"))
	expect_identical(test_obj@assay.used, "integrated")
	expect_identical(
		test_obj@call.string,
		"combineExpression(mini_contig_list, scRep_example)"
	)
	expect_equal(
		test_obj@params,
		list(
			cloneCall = "CTstrict", chain = "both", group.by = "none",  
			proportion = TRUE, filterNA = FALSE,
			cloneTypes = c(
				`None ( < X <= 0)` = 0, `Rare (0 < X <= 1e-04)` = 1e-04,
				`Small (1e-04 < X <= 0.001)` = 0.001,
				`Medium (0.001 < X <= 0.01)` = 0.01,
				`Large (0.01 < X <= 0.1)` = 0.1, 
				`Hyperexpanded (0.1 < X <= 1)` = 1
			),
			addLabel = FALSE
		)
	)
})
