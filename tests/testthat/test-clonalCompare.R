# test script for clonalCompare.R - testcases are NOT comprehensive!

test_that("clonalCompare works", {
  
  combined <- getCombined()

  expect_doppelganger(
    "clonalCompare_alluvial_plot",
    clonalCompare(
      combined, 
      top.clones = 10, 
      samples = c("P17B", "P17L"), 
      cloneCall="aa", 
      graph = "alluvial"
    )
  )
  
  expect_doppelganger(
    "clonalCompare_alluvial_order_plot",
    clonalCompare(
      combined, 
      top.clones = 10, 
      samples = c("P17B", "P17L"), 
      order.by = c("P17L", "P17B"),
      cloneCall="aa", 
      graph = "alluvial"
    )
  )
  
  expect_doppelganger(
    "clonalCompare_area_plot",
    clonalCompare(
      combined, 
      top.clones  = 10, 
      samples = c("P17B", "P17L"), 
      cloneCall="aa", 
      graph = "area")
  )
  
  expect_doppelganger(
  "clonalCompare_highlight.relabel_plot",
   clonalCompare(combined, 
                top.clones = 10,
                highlight.clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
                relabel.clones = TRUE,
                samples = c("P17B", "P17L"), 
                cloneCall="aa", 
                graph = "alluvial")
   )
  
})

test_that("clonalCompare works with exportTable and prop FALSE", {

  combined <- getCombined()

  getClonalCompareRes <- function(prop) {
    clonalCompare(
      combined,
      top.clones = 10,
      highlight.clones = c("CVVSDNTGGFKTIF_CASSVRRERANTGELFF", "NA_CASSVRRERANTGELFF"),
      relabel.clones = TRUE,
      samples = c("P17B", "P17L"),
      cloneCall = "aa",
      graph = "alluvial",
      exportTable = TRUE,
      proportion = prop
    )
  }

  countCompareRes <- getClonalCompareRes(prop = FALSE)
  propCompareRes <- getClonalCompareRes(prop = TRUE)

  expect_identical(
    countCompareRes %>% dplyr::select(-Count),
    propCompareRes %>% dplyr::select(-Proportion)
  )

  fullJoined <- getClonalCompareRes(FALSE) %>%
    dplyr::full_join(
      getClonalCompareRes(TRUE)
    )

  expect_setequal(
    colnames(fullJoined),
    c("clones", "Count", "original.clones", "Proportion", "Sample")
  )

  countPropFactor <- fullJoined$Count / fullJoined$Proportion

  expect_identical(
    as.integer(fullJoined$Count - fullJoined$Proportion * countPropFactor),
    integer(nrow(fullJoined))
  )
})
