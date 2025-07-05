# test script for positionalEntropy.R - testcases are NOT comprehensive!

combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))


test_that("Function returns a ggplot object by default", {
  p <- positionalEntropy(combined)
  expect_s3_class(p, "ggplot")
})

test_that("Function returns a data.frame when exportTable is TRUE", {
  df_out <- positionalEntropy(combined, exportTable = TRUE)
  expect_s3_class(df_out, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_out)))
})

test_that("Grouping with 'group.by' works correctly", {
  df_grouped <- positionalEntropy(combined, group.by = "Type", exportTable = TRUE)
  
  # Check for correct group names
  expect_equal(sort(unique(df_grouped$group)), c("B", "L"))
  
  # Check calculations for group B
  group_b_data <- subset(df_grouped, group == "B")
  expect_equal(subset(group_b_data, Position == 1)$entropy, 0) 
  expect_equal(subset(group_b_data, Position == 20)$entropy, 0) 
})

test_that("'method' parameter changes", {
  df_shannon <- positionalEntropy(combined, 
                                  method = "shannon", 
                                  exportTable = TRUE)
  expect_s3_class(df_shannon, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_shannon)))
  
  df_inv_simpson <- positionalEntropy(combined, 
                                      method = "inv.simpson", 
                                      exportTable = TRUE)
  expect_s3_class(df_inv_simpson, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_inv_simpson)))
  
  df_gini_simpson <- positionalEntropy(combined, 
                                       method = "gini.simpson", 
                                       exportTable = TRUE)
  expect_s3_class(df_gini_simpson, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_gini_simpson)))

  df_hill0 <- positionalEntropy(combined, 
                                method = "hill0", 
                                exportTable = TRUE)
  expect_s3_class(df_hill0 , "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_hill0 )))
  
  df_hill1 <- positionalEntropy(combined, 
                                method = "hill1", 
                                exportTable = TRUE)
  expect_s3_class(df_hill1, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_hill1)))
  
  df_hill2 <- positionalEntropy(combined, 
                                method = "hill2", 
                                exportTable = TRUE)
  expect_s3_class(df_hill2, "data.frame")
  expect_true(all(c("entropy", "Position", "group") %in% colnames(df_hill2)))

})

test_that("Ordering with 'order.by' works correctly", {
  order_vec <- c("L", "B")
  df_ordered <- positionalEntropy(combined, 
                                  group.by = "Type", 
                                  order.by = order_vec, 
                                  exportTable = TRUE)
  
  expect_s3_class(df_ordered$group, "factor")
  expect_equal(levels(df_ordered$group), order_vec)
})

t