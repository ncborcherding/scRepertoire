# test script for utils.R - testcases are NOT comprehensive!

test_that("'%!in%' works", {
    v <- c(1, 2, 3, 4, 5)
    
    expect_true(0 %!in% v)
    expect_true(6 %!in% v)
    expect_false(3 %!in% v)
    expect_false(5 %!in% v)
    expect_true(1 %!in% NULL)
    expect_true(list(1) %!in% NA)
})

# TODO off.the.chain
# TODO checkBlanks
# TODO groupList

test_that("checkList works", {
    data("contig_list")
    expect_identical(.checkList(contig_list), contig_list)
    expect_identical(.checkList(contig_list[[1]]), list(contig_list[[1]]))
    expect_identical(.checkList(contig_list[[1]])[[1]], contig_list[[1]])
    # no idea what to put to make the stop message happen. 
})

# TODO bound.input.return
# TODO get.coord
# TODO checkSingleObject
# TODO grabMeta

# Test cases
test_that("Modifying barcodes without ID works correctly", {
    samples <- c("sample1", "sample2")
    modified_data <- .modifyBarcodes(
        df = getdata("utils", "df_list"), samples = samples, ID = NULL
    )
    
    expected_modified_data <- list(
        data.frame(
            barcode = c("sample1_A", "sample1_B", "sample1_C"),
            value = c(10, 20, 30)
        ),
        data.frame(
            barcode = c("sample2_X", "sample2_Y", "sample2_Z"),
            value = c(100, 200, 300)
        )
    )
    
    expect_identical(modified_data, expected_modified_data)
})

##############
## resolved TODO: getdata is from tests/testthat/helper-testing_functions.R
##TODO annotate testthat functions???
test_that("Modifying barcodes with ID works correctly", {
    samples <- c("sample3", "sample4")
    ID <- c("id1", "id2")
    modified_data <- .modifyBarcodes(
        df = getdata("utils", "df_list"), samples = samples, ID = ID
    )
    
    expected_modified_data <- list(
        data.frame(
            barcode = c("sample3_id1_A", "sample3_id1_B", "sample3_id1_C"),
            value = c(10, 20, 30)
        ),
        data.frame(
            barcode = c("sample4_id2_X", "sample4_id2_Y", "sample4_id2_Z"),
            value = c(100, 200, 300)
        )
    )
    
    expect_identical(modified_data, expected_modified_data)
})


# TODO parseContigs
# TODO quiet

.quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}



test_that(".theCall works", {
  expect_equal(.theCall(NULL, "aa", check.df = FALSE), "CTaa")
  expect_equal(.theCall(NULL, "nt", check.df = FALSE), "CTnt")
  expect_equal(.theCall(NULL, "genes", check.df = FALSE), "CTgene")
  expect_equal(.theCall(NULL, "strict", check.df = FALSE), "CTstrict")
})
#TODO .theCall Add custom header

# TODO .constructConDfAndParseTCR !!!! Need to use testing data from the old version

test_that(".constructConDfAndParseTCR works", {
    # TODO create testdata with the original .parseTCR and test here. also do edgecases
})

# TODO .parseBCR
# TODO lengthDF
# TODO assignCT

test_that("makeGenes works for cellType T", {
    expect_identical(
        .makeGenes("T", getdata("utils", "makeGenes_T_input")),
        getdata("utils", "makeGenes_T_expected")
    )
})
# TODO makesGenes (cellType B)

# TODO short.check

# TODO select.gene

test_that("is_df_or_list_of_df works", {
    df <- data.frame(x = 1:5, y = letters[1:5])
    list_of_dfs <- list(
        data.frame(a = 1:3, b = letters[1:3]),
        data.frame(x = 4:6, y = letters[4:6])
    )
    mixed_list <- list(
        data.frame(a = 1:3, b = letters[1:3]),
        "not a dataframe"
    )
    
    expect_true(is_df_or_list_of_df(df))
    expect_true(is_df_or_list_of_df(list_of_dfs))
    expect_false(is_df_or_list_of_df(mixed_list))
    expect_false(is_df_or_list_of_df(list()))
    expect_false(is_df_or_list_of_df(c(1, 2, 3, 4, 5)))
})
