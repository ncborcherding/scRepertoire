# test script for percentGenes.R - testcases are NOT comprehensive!

test_that("percentGenes works", {

  combined <- getCombined()

  set.seed(42)
  expect_doppelganger(
    "percentGenes_plot",
    percentGenes(combined, 
                 chain = "TRB", 
                 gene = "V")
  )
  
  set.seed(42)
  expect_doppelganger(
    "percentGenes_order_plot",
    percentGenes(combined, 
                 chain = "TRB", 
                 gene = "V",
                 order.by = c("P17B","P18B","P19B","P20B","P17L","P18L","P19L","P20L"))
  )
  
  set.seed(42)
  expect_equal(
    percentGenes(combined, 
              chain = "TRB", 
              gene = "V",
              exportTable = TRUE),
      getdata("visualizations", "percentGenes_exportTable")
  )
})