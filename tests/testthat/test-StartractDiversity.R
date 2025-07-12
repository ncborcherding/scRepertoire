# test script for StartracDiversity.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)
scRep_example$Patient <- substring(scRep_example$orig.ident, 1, 3)
scRep_example$Type <- substring(scRep_example$orig.ident, 4, 4)

test_that("Output formats (table vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame
  table_output <- StartracDiversity(scRep_example, type = "Type", group.by = "Patient", exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Test default behavior returns a ggplot object
  plot_output <- StartracDiversity(scRep_example, type = "Type", group.by = "Patient")
  expect_s3_class(plot_output, "ggplot")
})

test_that("Data frame output has correct structure and content", {
  table_output <- StartracDiversity(scRep_example, type = "Type", group.by = "Patient", exportTable = TRUE)
  
  # Check for expected columns
  expected_cols <- c("group", "majorCluster", "migr", "tran", "expa")
  expect_true(all(expected_cols %in% names(table_output)))
  
  # Check number of groups
  num_patients <- length(unique(scRep_example$Patient))
  expect_equal(length(unique(table_output$group)), num_patients)
})

test_that("Plot output has correct structure", {
  plot_output <- StartracDiversity(scRep_example, type = "Type", group.by = "Patient")
  
  # Check for correct geometry and facetting
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomBoxplot")
  expect_s3_class(plot_output$facet, "FacetGrid")
  
  # Check labels
  expect_equal(plot_output$labels$y, "Index Score")
  expect_equal(plot_output$labels$x, "majorCluster")
})

test_that("`cloneCall` and `chain` parameters execute without error", {
  expect_silent(
    res_aa <- StartracDiversity(scRep_example, cloneCall = "aa", type = "Type", group.by = "Patient", exportTable = TRUE)
  )
  expect_s3_class(res_aa, "data.frame")
  expect_silent(
    res_tra <- StartracDiversity(scRep_example, chain = "TRA", type = "Type", group.by = "Patient", exportTable = TRUE)
  )
  expect_s3_class(res_tra, "data.frame")
})
