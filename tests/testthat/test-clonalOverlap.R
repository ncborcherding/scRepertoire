# test script for clonalOverlap.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list,
                       samples = c("P17B", "P17L", "P18B", "P18L",
                                   "P19B", "P19L", "P20B", "P20L"))
combined <- addVariable(combined, 
                        variable.name = "Type", 
                        variables = rep(c("B", "L"), 4))

test_that("Input validation for `method` works", {
  expect_error(clonalOverlap(combined, method = "invalid_method"), "Invalid method provided", fixed = TRUE)
})

test_that("Output formats (data.frame vs. plot) are correct", {
  # Test that exportTable = TRUE returns a data.frame
  df_output <- clonalOverlap(combined, method = "raw", exportTable = TRUE)
  expect_true(is.data.frame(df_output))
  
  # Test that the default behavior returns a ggplot object
  plot_output <- clonalOverlap(combined, method = "raw")
  expect_s3_class(plot_output, "ggplot")
})

methods <- c("raw", "jaccard", "overlap", "morisita", "cosine")
for (m in methods) {
    # Create a unique test name for each method
    test_that(paste("Method '", m, "' produces a valid matrix"), {
      mat <- clonalOverlap(combined, method = m, exportTable = TRUE)
      
      # Matrix should be square
      expect_equal(nrow(mat), ncol(mat))
      expect_equal(nrow(mat), length(combined))
      
      # Diagonal and lower triangle should be NA
      expect_true(all(is.na(mat[lower.tri(mat, diag = TRUE)])))
      
      # upper triangle should be numeric
      expect_true(is.numeric(mat[upper.tri(mat)]))
      
      # Check value ranges for indices
      if (m %in% c("jaccard", "overlap", "morisita", "cosine")) {
        expect_true(all(na.omit(mat) >= 0 & na.omit(mat) <= 1))
      } else if (m == "raw") {
        expect_true(all(na.omit(mat) >= 0))
      }
    })
}



methods <- c("raw", "jaccard", "morisita") 
for (m in methods) {
    test_that(paste("Grouping works for method '", m, "'"), {
      mat_grouped <- clonalOverlap(combined, group.by = "Type", method = m, exportTable = TRUE)
      
      # Matrix dimensions should match the number of groups
      expect_equal(nrow(mat_grouped), 2)
      expect_equal(ncol(mat_grouped), 2)
      
      # Check row and column names
      expect_equal(rownames(mat_grouped), c("B", "L"))
      expect_equal(colnames(mat_grouped), c("B", "L"))
    })
}


test_that("Calculations are correct (sanity check)", {
  # Create a simple, predictable list of data frames
  list_for_testing <- list(
    A = data.frame(CTaa = c("a", "b", "c")),
    B = data.frame(CTaa = c("b", "c", "d", "e")),
    C = data.frame(CTaa = c("x", "y"))
  )
  
  # Test raw count
  mat_raw <- clonalOverlap(list_for_testing, cloneCall = "aa", method = "raw", exportTable = TRUE)
  expect_equal(mat_raw["A", "B"], 2) # Overlap is "b", "c"
  expect_equal(mat_raw["A", "C"], 0) # No overlap
  
  # Test Jaccard index
  mat_jaccard <- clonalOverlap(list_for_testing, cloneCall = "aa", method = "jaccard", exportTable = TRUE)
  # Jaccard(A,B) = |A intersect B| / |A union B| = 2 / (3 + 4 - 2) = 2/5 = 0.4
  expect_equal(mat_jaccard["A", "B"], 0.4)
})

test_that("Plotting output has the correct components", {
  plot_output <- clonalOverlap(combined, method = "jaccard")
  
  # Check for correct geometry
  expect_s3_class(plot_output$layers[[1]]$geom, "GeomTile")
  
  # Check labels
  expect_equal(plot_output$labels$fill, "Jaccard")
})
