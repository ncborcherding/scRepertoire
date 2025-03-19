# test script for annotateInvariant.R - testcases are NOT comprehensive!
data("scRep_example")
test_obj <- combineExpression(getdata("combineContigs", "combined"), scRep_example)

test_that("annotateInvariant() handles correct input format", {
  # Mocking input data with required structure
  
  result <- annotateInvariant(test_obj, type = "MAIT", species = "human")
  
  expect_true(is.data.frame(result[[]]))
  expect_true(any(colnames(result[[]]) ==  "MAIT.score"))
})

test_that("annotateInvariant() fails with incorrect input format", {
  expect_error(annotateInvariant("wrong_input", type = "MAIT", species = "human"),
               "Please use the output of combineTCR() or combineExpression() as input.data")
})

test_that("annotateInvariant() handles species and type argument matching correctly", {
  expect_error(annotateInvariant(mock_input, type = "INVALID", species = "human"),
               "arg should be one of \"MAIT\", \"iNKT\"")
  
  expect_error(annotateInvariant(mock_input, type = "MAIT", species = "INVALID"),
               "arg should be one of \"mouse\", \"human\"")
})


test_that("annotateInvariant() returns zero scores for non-matching cells", {
  result <- annotateInvariant(test_obj, type = "MAIT", species = "human")
  
  expect_equal(sum(result$MAIT.score), 0) 
})

test_that("annotateInvariant() correctly integrates with Seurat object", {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    result <- annotateInvariant(test_obj, type = "MAIT", species = "human")
    
    expect_true("MAIT.score" %in% colnames(result@meta.data))
  }
})

