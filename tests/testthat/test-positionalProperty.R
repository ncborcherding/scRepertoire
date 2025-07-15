# test script for positionalProperty.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, 
                       samples = c("P17B", "P17L", "P18B", "P18L", 
                                   "P19B","P19L", "P20B", "P20L"))

test_that("positionalProperty: Output structure is correct", {
  
  # Test with exportTable = TRUE
  df_output <- positionalProperty(combined, 
                                  method = "kideraFactors", 
                                  exportTable = TRUE)
  
  expect_s3_class(df_output, "data.frame")
  
  expected_cols <- c("property", "position", "mean", "sd", "n", "ci_lower", "ci_upper", "group")
  expect_true(all(expected_cols %in% names(df_output)))
  
  # Test that it returns a ggplot object by default
  plot_output <- positionalProperty(combined, 
                                    method = "atchleyFactors")
  
  expect_s3_class(plot_output, "ggplot")
})

test_that("positionalProperty: Input validation works", {
  
  # Error on invalid method
  expect_error(
    positionalProperty(combined, 
                       method = "InvalidMethod"),
    "Please select a compatible method"
  )
  
})

test_that("positionalProperty: Core calculations are plausible", {
  
  df_output <- positionalProperty(combined, 
                                  chain = "TRB",
                                  aa.length = 15, 
                                  exportTable = TRUE)
  
  p17_b_subset <- df_output[df_output$group == "P17B", ]
  expect_true(all(p17_b_subset$n <= 2851)) 
})

test_that("positionalProperty: Parameters behave as expected", {
  
  # `aa.length` parameter
  df_len10 <- positionalProperty(combined, 
                                 aa.length = 10, 
                                 exportTable = TRUE)
  expect_equal(max(df_len10$position), 10)
  
  df_len25 <- positionalProperty(combined, 
                                 aa.length = 25, 
                                 exportTable = TRUE)
  expect_equal(max(df_len25$position), 25)
})

test_that("positionalProperty: ggplot object is correctly formed", {
  
  plot_output <- positionalProperty(combined, 
                                    method = "kideraFactors")
  
  # Check layers
  expect_true("GeomRibbon" %in% sapply(plot_output$layers, function(x) class(x$geom)[1]))
  expect_true("GeomLine" %in% sapply(plot_output$layers, function(x) class(x$geom)[2]))
  
  # Check aesthetics
  expect_equal(rlang::as_name(plot_output$mapping$y), "mean")
  expect_equal(rlang::as_name(plot_output$mapping$group), "group")
})
