combined <- combineTCR(
  contig_list,
  samples = c("P17B", "P17L", "P18B", "P18L", "P19B", "P19L", "P20B", "P20L")
)

test_that("percentKmer works for AAs", {
  top_30_aa_3mer_composition_matrix <- getdata(
    "percentKmer", "top_30_aa_3mer_composition_matrix"
  )

  expect_equal(
    percentKmer(combined, exportTable = TRUE),
    top_30_aa_3mer_composition_matrix
  )
})

test_that("percentKmer works for NTs", {
  top_30_nt_3mer_composition_matrix <- getdata(
    "percentKmer", "top_30_nt_3mer_composition_matrix"
  )

  expect_equal(
    percentKmer(combined, cloneCall = "nt", exportTable = TRUE),
    top_30_nt_3mer_composition_matrix
  )
})

# TODO test for cases where no kmers were counted (NA columns present)