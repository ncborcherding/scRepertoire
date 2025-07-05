# test script for percentAA.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))

test_that("Function returns a ggplot object by default", {
  p <- percentAA(combined)
  expect_s3_class(p, "ggplot")
})

test_that("Function returns a data.frame when exportTable is TRUE", {
  df_out <- percentAA(combined, exportTable = TRUE)
  expect_s3_class(df_out, "data.frame")
  # Check for expected columns
  expected_cols <- c("AminoAcid", "Position", "Frequency", "group")
  expect_true(all(expected_cols %in% colnames(df_out)))
})

test_that("Grouping with 'group.by' works correctly", {
  p_grouped <- percentAA(combined, group.by = "Type")
  df_grouped <- percentAA(combined, group.by = "Type", exportTable = TRUE)
  
  # Check data.frame for correct groups
  expect_equal(sort(unique(df_grouped$group)), c("B", "L"))
  
  # Check data integrity for group B
  group_b_data <- df_grouped[df_grouped$group == "B", ]
  pos1_C <- subset(group_b_data, Position == 1 & AminoAcid == "C")
  expect_equal(pos1_C$Frequency, 1) 
})


test_that("Ordering with 'order.by' works correctly", {
  # Groups are B, L. Default order is alphabetical. We want L, B.
  order_vec <- c("L", "B")
  df_ordered <- percentAA(combined, group.by = "Type", order.by = order_vec, exportTable = TRUE)
  
  # Check if the 'group' column is a factor with the specified levels
  expect_s3_class(df_ordered$group, "factor")
  expect_equal(levels(df_ordered$group), order_vec)
})
