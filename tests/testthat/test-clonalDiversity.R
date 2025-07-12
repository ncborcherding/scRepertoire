# test script for clonalDiversity.R - testcases are NOT comprehensive!

# Generate some data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                    "P19B","P19L", "P20B", "P20L"))

combined <- addVariable(combined, 
            variable.name = "Type", 
            variables = rep(c("B", "L"), 4))

test_that("clonalDiversity: Input validation works", {
  expect_error(
    clonalDiversity(combined_input_list, metric = "invalid_metric"),
    "Invalid `metric`"
  )
  expect_error(
    clonalDiversity(combined_input_list, metric = c("shannon", "inv.simpson")),
    "`metric` must be a single string."
  )
})

test_that("clonalDiversity: skip.boots = TRUE calculates correct values", {
  # Test with a list input
  result <- clonalDiversity(
    combined,
    cloneCall = "CTaa",
    metric = "shannon",
    skip.boots = TRUE,
    exportTable = TRUE
  )
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 8)
  expect_equal(result$value[result$Group == "P17B"], 4.108008, tolerance = 1e-6)
  expect_equal(result$value[result$Group == "P17L"], 7.33391, tolerance = 1e-6)
  
  # Test with a data.frame input
  result_df <- clonalDiversity(
    combined,
    cloneCall = "CTaa",
    metric = "inv.simpson",,
    skip.boots = TRUE,
    exportTable = TRUE
  )
  
  expect_equal(result_df$value[result_df$Group == "P17B"], 8.029998, tolerance = 1e-6)
  expect_equal(result_df$value[result_df$Group == "P17L"], 627.9127, tolerance = 1e-6)
})


test_that("clonalDiversity: Bootstrapping returns correct structure", {
  set.seed(42) # for reproducibility
  n_boots <- 10
  
  # Test with mean of bootstraps (default)
  result <- clonalDiversity(
    combined,
    cloneCall = "CTaa",
    metric = "shannon",
    n.boots = n_boots,
    skip.boots = FALSE,
    exportTable = TRUE
  )
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 8)
  expect_true(all(c("Group", "value", "metric") %in% names(result)))
  # Value should be a single mean value per group
  expect_true(is.numeric(result$value))
  
  # Test with return.boots = TRUE
  result_all_boots <- clonalDiversity(
    combined,
    cloneCall = "CTaa",
    metric = "shannon",
    n.boots = n_boots,
    return.boots = TRUE # This also sets exportTable = TRUE
  )
  expect_s3_class(result_all_boots, "data.frame")
  # Should have n.boots rows for each group
  expect_equal(nrow(result_all_boots), 8 * n_boots)
  expect_equal(sum(result_all_boots$Group == "P17B"), n_boots)
  expect_equal(sum(result_all_boots$Group == "P17L"), n_boots)
})


test_that("clonalDiversity: exportTable and plotting works", {
  # Returns a data.frame when exportTable = TRUE
  result_table <- clonalDiversity(
    combined,
    exportTable = TRUE
  )
  expect_s3_class(result_table, "data.frame")
  
  # Returns a ggplot object when exportTable = FALSE (default)
  result_plot <- clonalDiversity(
    combined,
    exportTable = FALSE
  )
  expect_s3_class(result_plot, "ggplot")
  
  # Test x.axis functionality returns a ggplot object
  result_plot_xaxis <- clonalDiversity(
    combined,
    x.axis = "Type",
    exportTable = FALSE
  )
  expect_s3_class(result_plot_xaxis, "ggplot")
})



