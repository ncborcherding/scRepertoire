
test_that("loadContigs works for BD, WAT3R, 10X, MiXCR, Immcantation, ParseBio, and Dandelion formats", {
  
  # BD format
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- loadContigs(BD, format = "BD")
  expected_bd <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_BD"))
  expect_identical(trial_bd, expected_bd)
  
  # WAT3R format
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- loadContigs(WAT3R, format = "WAT3R")
  expected_wat3r <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_WAT3R"))
  expect_identical(trial_wat3r, expected_wat3r)
  
  # 10X format (using pre-loaded contig_list data)
  data("contig_list")
  trial_10x <- loadContigs(contig_list[[1]], format = "10X")
  expected_10x <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_10x"))
  expect_identical(trial_10x, expected_10x)
  
  # MiXCR format
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- loadContigs(MIXCR, format = "MiXCR")
  expected_mixcr <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_MiXCR"))
  expect_identical(trial_mixcr, expected_mixcr)
  
  # Immcantation format
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- loadContigs(Immcantation, format = "Immcantation")
  expected_immcantation <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Immcantation"))
  expect_identical(trial_immcantation, expected_immcantation)
  
  # ParseBio format
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- loadContigs(Parse, format = "ParseBio")
  expected_parse <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Parse"))
  expect_identical(trial_parse, expected_parse)
  
  # Dandelion format
  Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
  trial_dandelion <- loadContigs(Dandelion, format = "Dandelion")
  expected_dandelion <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Dandelion"))
  expect_identical(trial_dandelion, expected_dandelion)
})

test_that("Auto-detection works for pre-loaded data", {
  # BD auto-detection
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd_auto <- loadContigs(BD, format = "auto")
  expected_bd <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_BD"))
  expect_identical(trial_bd_auto, expected_bd)
  
  # WAT3R auto-detection
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r_auto <- loadContigs(WAT3R, format = "auto")
  expected_wat3r <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_WAT3R"))
  expect_identical(trial_wat3r_auto, expected_wat3r)
  
  # 10X auto-detection (using pre-loaded contig_list)
  data("contig_list")
  trial_10x_auto <- loadContigs(contig_list[[1]], format = "auto")
  expected_10x <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_10x"))
  expect_identical(trial_10x_auto, expected_10x)
  
  # MiXCR auto-detection
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr_auto <- loadContigs(MIXCR, format = "auto")
  expected_mixcr <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_MiXCR"))
  expect_identical(trial_mixcr_auto, expected_mixcr)
  
  # Immcantation auto-detection
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation_auto <- loadContigs(Immcantation, format = "auto")
  expected_immcantation <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Immcantation"))
  expect_identical(trial_immcantation_auto, expected_immcantation)
  
  # ParseBio auto-detection
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse_auto <- loadContigs(Parse, format = "auto")
  expected_parse <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Parse"))
  expect_identical(trial_parse_auto, expected_parse)
  
  # Dandelion auto-detection
  Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
  trial_dandelion_auto <- loadContigs(Dandelion, format = "auto")
  expected_dandelion <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_Dandelion"))
  expect_identical(trial_dandelion_auto, expected_dandelion)
})

test_that("loadContigs works for TRUST4 format", {
  # TRUST4 file input test
  TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
  trial_trust4 <- loadContigs(TRUST4, format = "TRUST4")
  expected_trust4 <- .rmAllNaRowsFromLoadContigs(getdata("load", "loadContigs_TRUST4"))
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
  
  # Explicit format test
  result_airr <- loadContigs(tmp_dir_airr, format = "AIRR")
  expected_airr <- list(
    data.frame(
      barcode = "cellA",
      chain = "TRB",
      reads = as.integer(5),
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
  
  # Auto-detection test for AIRR directory
  result_airr_auto <- loadContigs(tmp_dir_airr, format = "auto")
  expect_identical(result_airr_auto, expected_airr)
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
