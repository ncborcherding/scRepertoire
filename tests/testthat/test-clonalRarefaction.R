# test script for clonalRarefaction.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))

combined <- combined[1:4] #subset to reduce size

test_that("Output formats (data.frame vs. plot) are correct", {
  table_output_list <- clonalRarefaction(combined, exportTable = TRUE, n.boots = 2)
  expect_type(table_output_list, "list")
  expect_s3_class(table_output_list$data, "data.frame")
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalRarefaction(combined, n.boots = 2)
  expect_s3_class(plot_output, "ggplot")
})

test_that("`plot.type` parameter correctly generates plots", {
  expect_silent(p1 <- clonalRarefaction(combined, plot.type = 1, n.boots = 2))
  expect_s3_class(p1, "ggplot")
  expect_silent(p2 <- clonalRarefaction(combined, plot.type = 2, n.boots = 2))
  expect_s3_class(p2, "ggplot")
  expect_silent(p3 <- clonalRarefaction(combined, plot.type = 3, n.boots = 2))
  expect_s3_class(p3, "ggplot")
})

test_that("`hill.numbers` parameter is passed to iNEXT correctly", {
  df_q0 <- clonalRarefaction(combined, hill.numbers = 0, n.boots = 2, exportTable = TRUE)$data
  expect_true(all(df_q0$Order.q == 0))
  df_q01 <- clonalRarefaction(combined, hill.numbers = c(0, 1), n.boots = 2, exportTable = TRUE)$data
  expect_true(all(unique(df_q01$Order.q) %in% c(0, 1)))
})

test_that("`group.by` parameter correctly aggregates data", {
  # Group by the "Type" variable
  df_grouped <- clonalRarefaction(combined, group.by = "Type", n.boots = 2, exportTable = TRUE)$data
  num_groups <- 2
  expect_equal(length(unique(df_grouped$Assemblage)), num_groups)
  expect_true(all(unique(df_grouped$Assemblage) %in% c("B", "L")))
})

test_that("`cloneCall` and `chain` parameters execute correctly", {
  expect_silent(
    res_aa <- clonalRarefaction(combined, cloneCall = "aa", n.boots = 2, exportTable = TRUE)
  )
  expect_type(res_aa, "list")
  expect_silent(
    res_tra <- clonalRarefaction(combined, chain = "TRA", n.boots = 2, exportTable = TRUE)
  )
  expect_type(res_tra, "list")
})


