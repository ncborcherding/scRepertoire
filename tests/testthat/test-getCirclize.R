# test script for getCirclize.R - testcases are NOT comprehensive!

# Data to use
combined <- combineTCR(contig_list, samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L"))
scRep_example <- combineExpression(combined, scRep_example)

test_that("Output has the correct format", {
  df <- getCirclize(scRep_example, group.by = "seurat_clusters")
  
  expect_s3_class(df, "data.frame")
  expect_equal(names(df), c("from", "to", "value"))
  expect_true(is.numeric(df$value))
})

test_that("`include.self = FALSE` correctly removes self-referential rows", {
  df <- getCirclize(scRep_example, group.by = "seurat_clusters", include.self = FALSE)
  self_referential_rows <- df[df$from == df$to, ]
  expect_equal(nrow(self_referential_rows), 0)
  expect_gt(nrow(df), 0)
})


test_that("Function handles different `cloneCall` methods", {
  expect_silent(
    df_strict <- getCirclize(scRep_example, group.by = "seurat_clusters", cloneCall = "strict")
  )
  df_aa <- getCirclize(scRep_example, group.by = "seurat_clusters", cloneCall = "aa")
  expect_equal(df_strict$value, df_aa$value)
})

