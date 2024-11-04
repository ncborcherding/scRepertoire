# test script for clonalCluster.R - testcases are NOT comprehensive!

getTestCombinedSeuratObj <- function() {
  test_obj <- combineExpression(getCombined(), get(data("scRep_example")))
  test_obj$Patient <- substr(test_obj$orig.ident,1,3)
  test_obj$Type <- substr(test_obj$orig.ident,4,4)
  test_obj
}

test_that("clonalCluster works", {
 
  test_obj <- getTestCombinedSeuratObj()
  combined <- getCombined()
  withr::local_seed(42)

  expect_equal(
    clonalCluster(combined[[1]], 
                  chain = "TRB", 
                  sequence = "aa"),
    getdata("clustering", "clonalCluster_TRBaa_data")
  )
  
  set.seed(42)
  expect_equal(
  clonalCluster(combined[1:2], 
                chain = "TRB", 
                sequence = "aa"),
  getdata("clustering", "clonalCluster_2sample_data")
  )
  
  set.seed(42)
  test_obj <- clonalCluster(test_obj, 
                      chain = "TRB", 
                      sequence = "aa", 
                      group.by = "Patient")
  expect_equal(
    test_obj@meta.data[, 7:16],
    getdata("clustering", "clonalCluster_TRBaa_metadata")
  )
  
  
  BCR <- read.csv("https://www.borch.dev/uploads/contigs/b_contigs.csv")
  combined.BCR <- combineBCR(BCR, 
                             samples = "S1")
  expect_equal(
      clonalCluster(combined.BCR, 
                    chain = "IGH", 
                    sequence = "aa"),
      getdata("clustering", "clonalCluster_IGHaa_data")
  )
  
})

# test_that("clonalCluster works with custom threshold", { # fails atm

#   test_obj <- getTestCombinedSeuratObj()
#   withr::local_seed(42)

#   colorblind_vector <- hcl.colors(n=7, palette = "inferno", fixup = TRUE)

# test_obj <- clonalCluster(test_obj, 
#                                chain = "TRA", 
#                                sequence = "aa", 
#                                threshold = 0.85, 
#                                group.by = "Patient")

# Seurat::DimPlot(scRep_example, group.by = "TRA_cluster") +
#     scale_color_manual(values =  hcl.colors(n=length(unique(scRep_example@meta.data[,"TRA_cluster"])), "inferno")) + 
#   Seurat::NoLegend() + 
#   theme(plot.title = element_blank()) %>%
#   print()

# })

#TODO Add exportgraph test