test_that("percentKmer works for AAs", {
  top_30_aa_3mer_composition_matrix <- getdata(
    "percentKmer", "top_30_aa_3mer_composition_matrix"
  )

  expect_equal(
    percentKmer(getCombined(), cloneCall = "aa", exportTable = TRUE),
    top_30_aa_3mer_composition_matrix
  )
})

test_that("percentKmer works for NTs", {
  top_30_nt_3mer_composition_matrix <- getdata(
    "percentKmer", "top_30_nt_3mer_composition_matrix"
  )

  expect_equal(
    percentKmer(getCombined(), cloneCall = "nt", exportTable = TRUE),
    top_30_nt_3mer_composition_matrix
  )
  
  combined <- getCombined()
  combined <- addVariable(combined, 
                          variable.name = "Type", 
                          variables = rep(c("B", "L"), 4))
  expect_doppelganger(
    "percentKmer_group_motif2_plot",
    percentKmer(combined, 
                  motif.length = 2,
                  cloneCall = "aa",
                  group.by = "Type")
  )
})

# TODO test for cases where no kmers were counted (NA columns present)