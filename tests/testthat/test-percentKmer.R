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

test_that("tokenize_sequence works", {
  expect_equal(.tokenize_sequence("CAYRSAQAGGTSYGKLTF", 3),
               c("CAY", "AYR", "YRS", "RSA", "SAQ","AQA", "QAG", 
                 "AGG", "GGT", "GTS", "TSY", "SYG", "YGK", "GKL", "KLT", "LTF")
  )
})

test_that("tokenize_multiple_sequences works", {
  expect_equal(.tokenize_multiple_sequences(c("TESTING", "MULTIPLE", "SEQUENCE", "TOKENIZER"), 4),
               list(TESTING = c("TEST", "ESTI", "STIN", "TING"),
                    MULTIPLE = c("MULT", "ULTI", "LTIP", "TIPL", "IPLE"),
                    SEQUENCE = c("SEQU", "EQUE", "QUEN", "UENC", "ENCE"),
                    TOKENIZER = c("TOKE", "OKEN", "KENI", "ENIZ", "NIZE", "IZER"))
  )
                    
})

# TODO test for cases where no kmers were counted (NA columns present)