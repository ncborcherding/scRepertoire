# test script for quietVDJgenes.R - testcases are NOT comprehensive!

# Sample gene vectors for testing
tcr_genes <- c("TRAV1-1", "TRBV2", "TRGV9", "TRDV3", "TRBJ1-1", "TRDJ2")
bcr_genes <- c("IGHV3-23", "IGKV1-5", "IGLV3-1", "IGHJ4", "IGLC2", "JCHAIN")
other_genes <- c("CD3D", "GZMB", "MS4A1", "CD19")
pseudogenes <- c("IGHV1-12", "IGKV1-22") # From getHumanIgPseudoGenes()
all_genes <- c(tcr_genes, bcr_genes, other_genes, pseudogenes)


test_that("shouldQuietTcrGene correctly identifies TCR genes", {
  expected <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  expect_equal(shouldQuietTcrGene(tcr_genes), expected)
  
  expect_false(any(shouldQuietTcrGene(other_genes)))
  expect_false(any(shouldQuietTcrGene(bcr_genes)))
})

test_that("shouldQuietBcrGene correctly identifies BCR genes", {
  # All BCR genes except pseudogenes, plus JCHAIN, should be quieted
  bcr_genes_to_quiet <- c("IGHV3-23", "IGKV1-5", "IGLV3-1", "IGHJ4", "IGLC2", "JCHAIN")
  expect_true(all(shouldQuietBcrGene(bcr_genes_to_quiet)))
  expect_false(any(shouldQuietBcrGene(pseudogenes)))
  expect_false(any(shouldQuietBcrGene(other_genes)))
})

test_that("getHumanIgPseudoGenes returns expected character vector", {
  pseudos <- getHumanIgPseudoGenes()
  expect_type(pseudos, "character")
  expect_gt(length(pseudos), 0)
  expect_true("IGHV1-12" %in% pseudos)
  expect_true("IGLV1-41" %in% pseudos)
})

test_that("quietTCRgenes.default removes only TCR genes", {
  cleaned_genes <- quietTCRgenes.default(all_genes)
  expect_length(cleaned_genes, length(bcr_genes) + length(other_genes) + length(pseudogenes))
  expect_false(any(tcr_genes %in% cleaned_genes))
  expect_true(all(other_genes %in% cleaned_genes))
  expect_true(all(bcr_genes %in% cleaned_genes))
})

test_that("quietBCRgenes.default removes only BCR genes and JCHAIN", {
  cleaned_genes <- quietBCRgenes.default(all_genes)
  expected_genes <- c(tcr_genes, other_genes, pseudogenes)
  expect_equal(sort(cleaned_genes), sort(expected_genes))
  expect_false("JCHAIN" %in% cleaned_genes)
})

test_that("quietVDJgenes works as a wrapper on character vectors", {
  cleaned_genes <- quietVDJgenes(all_genes)
  expected_genes <- c(other_genes, pseudogenes)
  expect_equal(sort(cleaned_genes), sort(expected_genes))
  expect_length(cleaned_genes, length(other_genes) + length(pseudogenes))
})

test_that("Default methods handle empty input", {
  expect_equal(quietTCRgenes.default(character(0)), character(0))
  expect_equal(quietBCRgenes.default(character(0)), character(0))
  expect_equal(quietVDJgenes(character(0)), character(0))
})



