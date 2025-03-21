library(testthat)
library(yourPackageName)  # replace with your package name if needed

context("Testing loadContigs function")

test_that("loadContigs works for BD, WAT3R, 10X, MiXCR, Immcantation, ParseBio, and Dandelion formats", {
  
  # BD format
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- loadContigs(BD, format = "BD")
  expected_bd <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_BD"))
  expect_identical(trial_bd, expected_bd)
  
  # WAT3R format
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- loadContigs(WAT3R, format = "WAT3R")
  expected_wat3r <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_WAT3R"))
  expect_identical(trial_wat3r, expected_wat3r)
  
  # 10X format (using pre-loaded contig_list data)
  data("contig_list")
  trial_10x <- loadContigs(contig_list[[1]], format = "10X")
  expected_10x <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_10x"))
  expect_identical(trial_10x, expected_10x)
  
  # MiXCR format
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- loadContigs(MIXCR, format = "MiXCR")
  expected_mixcr <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_MiXCR"))
  expect_identical(trial_mixcr, expected_mixcr)
  
  # Immcantation format
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- loadContigs(Immcantation, format = "Immcantation")
  expected_immcantation <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Immcantation"))
  expect_identical(trial_immcantation, expected_immcantation)
  
  # ParseBio format
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- loadContigs(Parse, format = "ParseBio")
  expected_parse <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Parse"))
  expect_identical(trial_parse, expected_parse)
  
  # Dandelion format
  Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
  trial_dandelion <- loadContigs(Dandelion, format = "Dandelion")
  expected_dandelion <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Dandelion"))
  expect_identical(trial_dandelion, expected_dandelion)
})

test_that("loadContigs works for TRUST4 format", {
  # TRUST4 file input test
  TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
  trial_trust4 <- loadContigs(TRUST4, format = "TRUST4")
  expected_trust4 <- rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_TRUST4"))
  expect_identical(trial_trust4, expected_trust4)
  
  # One-row TRUST4 test using a minimal data frame
  oneRowTrust4Input <- data.frame(
    "#barcode" = "CGTAGCGGTGATAAGT-1",
    cell_type = "B",
    chain1 = "*",
    chain2 = "IGKV1D-43,*,IGKJ1,IGKC,TGTCAACAGTATAGTAGTGTCCCCTGGACGTTC,CQQYSSVPWTF,6.00,CGTAGCGGTGATAAGT-1_2,76.00,0",
    secondary_chain1 = "*",
    secondary_chain2 = "*",
    stringsAsFactors = FALSE
  )
  expected_parsed_trust4 <- list(
    data.frame(
      barcode = "CGTAGCGGTGATAAGT-1",
      v_gene = "IGKV1D-43",
      d_gene = "None",
      j_gene = "IGKJ1",
      c_gene = "IGKC",
      cdr3_nt = "TGTCAACAGTATAGTAGTGTCCCCTGGACGTTC",
      cdr3 = "CQQYSSVPWTF",
      reads = "6.00",
      chain = "IGK",
      stringsAsFactors = FALSE
    )
  )
  expect_identical(
    loadContigs(oneRowTrust4Input, format = "TRUST4"),
    expected_parsed_trust4
  )
})

test_that("loadContigs works with JSON input (directory mode)", {
  # Create a temporary JSON file with minimal data
  tmp_json <- tempfile(fileext = ".json")
  json_content <- '[{"cell_id": "cell1", "locus": "TRA", "consensus_count": 10, "v_call": "TRAV1", "d_call": "TRAD1", "j_call": "TRAJ1", "c_call": "TRAC1", "junction": "ATGC", "junction_aa": "M"}]'
  writeLines(json_content, tmp_json)
  
  # Copy to a temporary directory
  tmp_dir <- file.path(tempdir(), "json_test")
  dir.create(tmp_dir, showWarnings = FALSE)
  file.copy(tmp_json, file.path(tmp_dir, "test_data.json"), overwrite = TRUE)
  
  result_json <- loadContigs(tmp_dir, format = "JSON")
  expected_json <- list(
    data.frame(
      barcode = "cell1",
      chain = "TRA",
      reads = 10,
      v_gene = "TRAV1",
      d_gene = "TRAD1",
      j_gene = "TRAJ1",
      c_gene = "TRAC1",
      cdr3_nt = "ATGC",
      cdr3 = "M",
      stringsAsFactors = FALSE
    )
  )
  expect_identical(result_json, expected_json)
})

test_that("loadContigs works with AIRR input (directory mode)", {
  # Create a temporary TSV file for AIRR format
  tmp_tsv <- tempfile(fileext = ".tsv")
  airr_data <- data.frame(
    cell_id = "cellA",
    locus = "TRB",
    consensus_count = 5,
    v_call = "TRBV1",
    d_call = "TRBD1",
    j_call = "TRBJ1",
    c_call = "TRBC1",
    junction = "ATGCGT",
    junction_aa = "ML",
    stringsAsFactors = FALSE
  )
  write.table(airr_data, tmp_tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Copy to a temporary directory with proper filename
  tmp_dir_airr <- file.path(tempdir(), "airr_test")
  dir.create(tmp_dir_airr, showWarnings = FALSE)
  file.copy(tmp_tsv, file.path(tmp_dir_airr, "airr_rearrangement.tsv"), overwrite = TRUE)
  
  result_airr <- loadContigs(tmp_dir_airr, format = "AIRR")
  expected_airr <- list(
    data.frame(
      barcode = "cellA",
      chain = "TRB",
      reads = 5,
      v_gene = "TRBV1",
      d_gene = "TRBD1",
      j_gene = "TRBJ1",
      c_gene = "TRBC1",
      cdr3_nt = "ATGCGT",
      cdr3 = "ML",
      stringsAsFactors = FALSE
    )
  )
  expect_identical(result_airr, expected_airr)
})

test_that("loadContigs works with directory input for 10X format", {
  # Create a temporary CSV file mimicking a 10X contig file
  tmp_csv <- tempfile(fileext = ".csv")
  df_10x <- data.frame(
    barcode = c("cell1", "cell2"),
    chain = c("TRA", "TRB"),
    productive = c(TRUE, TRUE),
    cdr3 = c("ATGC", "TACG"),
    reads = c(100, 200),
    stringsAsFactors = FALSE
  )
  write.csv(df_10x, tmp_csv, row.names = FALSE)
  
  # Copy the file into a temporary directory with the expected name
  tmp_dir_10x <- file.path(tempdir(), "tenX_test")
  dir.create(tmp_dir_10x, showWarnings = FALSE)
  target_file <- file.path(tmp_dir_10x, "filtered_contig_annotations.csv")
  file.copy(tmp_csv, target_file, overwrite = TRUE)
  
  # Load from directory; note that filtering in .parse10x will remove non-productive chains and chains with "None"
  # For this test, we expect our minimal data to be returned (after ordering/sanitization)
  result_10x <- loadContigs(tmp_dir_10x, format = "10X")
  expected_10x_dir <- list(order_df(sanitize_empty(df_10x)))
  expect_identical(result_10x, expected_10x_dir)
})

test_that("loadContigs errors on unsupported format and invalid input types", {
  # Unsupported format
  expect_error(loadContigs("dummy", format = "UnsupportedFormat"))
  expect_error(loadContigs(12345, format = "10X"))
})

test_that("loadContigs returns empty list when no matching files are found in a directory", {
  empty_dir <- file.path(tempdir(), "empty_test_dir")
  dir.create(empty_dir, showWarnings = FALSE)
  # Ensure the directory is empty
  file.remove(list.files(empty_dir, full.names = TRUE))
  expect_warning(res <- loadContigs(empty_dir, format = "10X"))
  expect_identical(res, list())
})
