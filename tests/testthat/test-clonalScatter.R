# test script for clonalScatter.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined,
                        variable.name = "Type",
                        variables = rep(c("B", "L"), 4))

test_that("Output formats (table vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame
  table_output <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", exportTable = TRUE)
  expect_s3_class(table_output, "data.frame")
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L")
  expect_s3_class(plot_output, "ggplot")
})

test_that("`graph` parameter correctly sets up data and plot", {
  plot_prop <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", graph = "proportion")
  table_prop <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", graph = "proportion", exportTable = TRUE)
  expect_true(all(c("P17B.fraction", "P17L.fraction") %in% names(table_prop)))
  expect_true(any(sapply(plot_prop$layers, function(x) inherits(x$geom, "GeomAbline"))))
  expect_s3_class(plot_prop$scales$scales[[2]], "ScaleContinuous") # x-axis scale
  expect_s3_class(plot_prop$scales$scales[[3]], "ScaleContinuous") # y-axis scale
  expect_equal(plot_prop$scales$scales[[2]]$trans$name, "sqrt")
  
  # Test with graph = "count"
  plot_count <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", graph = "count")
  table_count <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", graph = "count", exportTable = TRUE)
  expect_true(all(c("P17B.fraction", "P17L.fraction") %in% names(table_count)))
  expect_false(any(sapply(plot_count$layers, function(x) inherits(x$geom, "GeomAbline"))))
})

test_that("Clone classification logic is correct", {
  list_for_testing <- list(
    X1 = data.frame(CTaa = c("A", "B", "C", "C", "C")), # A=1, B=1, C=3
    X2 = data.frame(CTaa = c("B", "D", "D")),          # B=1, D=2
    X3 = data.frame(CTaa = c("E"))
  )
  
  df <- clonalScatter(list_for_testing, cloneCall = "aa", x.axis = "X1", y.axis = "X2", exportTable = TRUE)
  clone_A <- df[df$Var1 == "A", ]
  clone_B <- df[df$Var1 == "B", ]
  clone_C <- df[df$Var1 == "C", ]
  clone_D <- df[df$Var1 == "D", ]
  expect_equal(clone_A$class, "X1.singlet")   # Count of 1 in X1, 0 in X2
  expect_equal(clone_C$class, "X1.expanded")  # Count > 1 in X1, 0 in X2
  expect_equal(clone_D$class, "X2.expanded")  # Count > 1 in X2, 0 in X1
  expect_equal(clone_B$class, "dual.expanded") # Present in both
})

test_that("`dot.size` parameter works correctly", {
  # Test with default dot.size = "total"
  df_total <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", dot.size = "total", exportTable = TRUE)
  expect_true("size" %in% names(df_total))
  expect_equal(df_total$size, df_total$P17B + df_total$P17L)
  
  # Test with a third sample for dot size
  df_third <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L", dot.size = "P18B", exportTable = TRUE)
  expect_true("P18B" %in% names(df_third))
})

test_that("`group.by` parameter correctly aggregates data", {
  # Group by the "Type" variable
  df_grouped <- clonalScatter(combined, cloneCall="aa", group.by = "Type", x.axis = "B", y.axis = "L", exportTable = TRUE)
  
  # Manually aggregate data for comparison
  all_data <- do.call(rbind, combined)
  b_clones <- subset(all_data, Type == "B")
  l_clones <- subset(all_data, Type == "L")
  b_counts <- as.data.frame(table(b_clones$CTaa))
  l_counts <- as.data.frame(table(l_clones$CTaa))
  shared_clone <- intersect(b_counts$Var1, l_counts$Var1)[1]
  manual_b_count <- b_counts[b_counts$Var1 == shared_clone, "Freq"]
  manual_l_count <- l_counts[l_counts$Var1 == shared_clone, "Freq"]
  func_counts <- df_grouped[df_grouped$Var1 == shared_clone, ]
  expect_equal(func_counts$B, manual_b_count)
  expect_equal(func_counts$L, manual_l_count)
})

test_that("Plotting output has the correct components", {
  plot_output <- clonalScatter(combined, x.axis = "P17B", y.axis = "P17L")
  
  # Check for correct geometry
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomAbline") # First layer is abline for proportion
  expect_s3_class(plot_output$layers[[2]]$geom, "GeomPoint")
  
  # Check labels
  expect_equal(plot_output$labels$x, "P17B")
  expect_equal(plot_output$labels$y, "P17L")
  expect_equal(plot_output$labels$fill, "class")
  expect_equal(plot_output$labels$size, "Total n")
})
