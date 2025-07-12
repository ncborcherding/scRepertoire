# test script for clonalQuant.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Input validation for `group.by` works", {
  error_message <- "Only one item in the group.by variable can be listed."
  expect_error(clonalQuant(combined, group.by = c("Type", "Patient")),
               error_message, fixed = TRUE)
})

test_that("Output formats (table vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame
  table_output <- clonalQuant(combined, exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalQuant(combined)
  expect_s3_class(plot_output, "ggplot")
})

test_that("`scale` parameter correctly calculates values", {
  table_raw <- clonalQuant(combined, exportTable = TRUE, scale = FALSE)
  expect_true(all(c("contigs", "values", "total") %in% names(table_raw)))
  expect_true(is.numeric(table_raw$contigs))
  expect_false("scaled" %in% names(table_raw))
  table_scaled <- clonalQuant(combined, exportTable = TRUE, scale = TRUE)
  expect_true("scaled" %in% names(table_scaled))
  expect_equal(table_scaled$scaled, (table_scaled$contigs / table_scaled$total) * 100)
  expect_true(all(table_scaled$scaled >= 0 & table_scaled$scaled <= 100))
})

test_that("Data frame contents are correct", {
  table_output <- clonalQuant(combined, cloneCall = "aa", exportTable = TRUE)
  expect_equal(nrow(table_output), length(combined))
  first_sample_data <- combined[[1]]
  manual_unique_count <- length(na.omit(unique(first_sample_data$CTaa)))
  manual_total_count <- length(na.omit(first_sample_data$CTaa))
  expect_equal(table_output$contigs[1], manual_unique_count)
  expect_equal(table_output$total[1], manual_total_count)
})

test_that("`group.by` parameter functions correctly", {
  table_grouped <- clonalQuant(combined, group.by = "Type", exportTable = TRUE)
  expect_equal(nrow(table_grouped), 2) # Should have 2 groups: "B" and "L"
  expect_true("Type" %in% names(table_grouped))
  grouped_list <- scRepertoire:::.groupList(combined, "Type")
  b_group_data <- grouped_list$B
  manual_unique_b <- length(na.omit(unique(b_group_data$CTstrict)))
  manual_total_b <- length(na.omit(b_group_data$CTstrict))
  b_row <- table_grouped[table_grouped$Type == "B",]
  expect_equal(b_row$contigs, manual_unique_b)
  expect_equal(b_row$total, manual_total_b)
})

test_that("Plotting output has the correct components and labels", {
  plot_raw <- clonalQuant(combined, scale = FALSE)
  expect_s3_class(plot_raw, "ggplot")
  expect_equal(plot_raw$labels$y, "Unique Clones")
  expect_equal(plot_raw$labels$fill, "Samples")
  plot_scaled <- clonalQuant(combined, scale = TRUE)
  expect_s3_class(plot_scaled, "ggplot")
  expect_equal(plot_scaled$labels$y, "Percent of Unique Clones")
  plot_grouped <- clonalQuant(combined, group.by = "Type")
  expect_s3_class(plot_grouped, "ggplot")
  expect_equal(plot_grouped$labels$fill, "Type")
})
