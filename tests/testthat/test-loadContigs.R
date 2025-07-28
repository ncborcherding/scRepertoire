# test script for loadContigs.R - testcases are NOT comprehensive!

check_loadContigs_output <- function(loaded_data) {
  # Check if the output is a list containing a single data frame
  expect_type(loaded_data, "list")
  expect_length(loaded_data, 1)
  df <- loaded_data[[1]]
  expect_s3_class(df, "data.frame")
  
  # Check if the data frame is not empty
  expect_gt(nrow(df), 0)
  
  # Check for the presence of essential standardized columns
  expected_cols <- c("barcode", "chain", "reads", "v_gene", 
                     "d_gene", "j_gene", "c_gene", "cdr3_nt", "cdr3")
  expect_true(all(expected_cols %in% names(df)))
  
  # Check data types of key columns
  expect_type(df$chain, "character")
  expect_type(df$cdr3, "character")
  # Reads should be numeric/integer after parsing
  expect_true(is.numeric(df$reads) || is.integer(df$reads))
}


test_that("loadContigs correctly processes various formats from URL", {
  #TRUST4 format
  TRUST4 <- read.csv("https://www.borch.dev/uploads/contigs/TRUST4_contigs.csv")
  trial_trust4 <- loadContigs(TRUST4, format = "TRUST4")
  check_loadContigs_output(trial_trust4)
  
  # BD format
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- loadContigs(BD, format = "BD")
  check_loadContigs_output(trial_bd)
  
  # WAT3R format
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- loadContigs(WAT3R, format = "WAT3R")
  check_loadContigs_output(trial_wat3r)
  
  # 10X format (using pre-loaded contig_list data)
  data("contig_list")
  trial_10x <- loadContigs(contig_list[[1]], format = "10X")
  check_loadContigs_output(trial_10x)
  
  # MiXCR format
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- loadContigs(MIXCR, format = "MiXCR")
  check_loadContigs_output(trial_mixcr)
  
  # Immcantation format
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- loadContigs(Immcantation, format = "Immcantation")
  check_loadContigs_output(trial_immcantation)
  
  # ParseBio format
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- loadContigs(Parse, format = "ParseBio")
  check_loadContigs_output(trial_parse)
  
  # Dandelion format
  Dandelion <- read.csv("https://www.borch.dev/uploads/contigs/Dandelion_contigs.csv")
  trial_dandelion <- loadContigs(Dandelion, format = "Dandelion")
  check_loadContigs_output(trial_dandelion)
})

test_that("loadContigs correctly auto-detects and processes various formats", {
  
  BD <- read.csv("https://www.borch.dev/uploads/contigs/BD_contigs.csv")
  trial_bd <- expect_message(
    loadContigs(BD, format = "auto"), 
    "Automatically detected format from data: AIRR"
  )
  
  WAT3R <- read.csv("https://www.borch.dev/uploads/contigs/WAT3R_contigs.csv")
  trial_wat3r <- expect_message(
    loadContigs(WAT3R, format = "auto"), 
    "Automatically detected format from data: WAT3R"
  )
  
  trial_10x <- expect_message(
      loadContigs(contig_list[[1]], format = "auto"), 
      "Automatically detected format from data: 10X"
  )
  
  MIXCR <- read.csv("https://www.borch.dev/uploads/contigs/MIXCR_contigs.csv")
  trial_mixcr <- expect_message(
    loadContigs(MIXCR, format = "auto"), 
    "Automatically detected format from data: MiXCR"
  )
  
  Immcantation <- read.csv("https://www.borch.dev/uploads/contigs/Immcantation_contigs.csv")
  trial_immcantation <- expect_message(
    loadContigs(Immcantation, format = "auto"), 
    "Automatically detected format from data: Immcantation"
  )
  
  Parse <- read.csv("https://www.borch.dev/uploads/contigs/Parse_contigs.csv")
  trial_parse <- expect_message(
    loadContigs(Parse, format = "auto"), 
    "Automatically detected format from data: ParseBio"
  )
})

test_that("loadContigs works with AIRR input (directory mode)", {
  # This test remains unchanged as it is self-contained and does not use internal data.
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
})

test_that("loadContigs returns empty list when no matching files are found in a directory", {
  empty_dir <- file.path(tempdir(), "empty_test_dir")
  dir.create(empty_dir, showWarnings = FALSE)
  file.remove(list.files(empty_dir, full.names = TRUE))
  expect_warning(res <- loadContigs(empty_dir, format = "10X"))
  expect_identical(res, list())
})
