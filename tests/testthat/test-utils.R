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
# TODO groupList
# TODO checkSingleObject
# TODO .parseBCR
# TODO lengthDF



test_that(".list.input.return works", {
    data("scRep_example")
    test_obj <- combineExpression(getCombined(), scRep_example)
    expect_equal(.list.input.return(test_obj, NULL),
                 getdata("utils", "list.input.return_data"))
    
})

test_that(".bound.input.return works", {
  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  expect_equal(.bound.input.return(test_obj),
               getdata("utils", "bound.input.return_data"))
  
})

test_that(".grabMeta works", {
  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  expect_equal(.grabMeta(test_obj),
               getdata("utils", "grabMeta_data"))
  
})

test_that(".expression2List works", {
  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  expect_equal(.expression2List(test_obj, NULL),
               getdata("utils", "expression2List_data"))
  expect_equal(.expression2List(test_obj, "orig.ident"),
               getdata("utils", "expression2List_by.orig.ident_data"))
  
})

test_that(".checkList works", {
  data("contig_list")
  expect_identical(.checkList(contig_list), contig_list)
  expect_identical(.checkList(contig_list[[1]]), list(contig_list[[1]]))
  expect_identical(.checkList(contig_list[[1]])[[1]], contig_list[[1]])
  # no idea what to put to make the stop message happen. 
})

test_that(".data.wrangle works", {
  data("scRep_example")
  test_obj <- combineExpression(getCombined(), scRep_example)
  expect_equal(.data.wrangle(test_obj, split.by = NULL, cloneCall = "CTaa", chain = "both"),
               getdata("utils", "data.wrangle_data"))
})


test_that(".checkContigs works", {
    data("contig_list")
    expect_equal(.checkContigs(contig_list), 
                 getdata("utils", "checkContigs_data"))
})

test_that(".checkBlanks works", {
  list_of_items.1 <- list(df1 = data.frame(x = 1:3, CTaa = letters[1:3]), 
                        df2 = data.frame(x = 6:9, CTaa = rep(NA,4)))
  expect_equal(.checkBlanks(list_of_items.1, "CTaa"), 
               list(df1 = data.frame(x = 1:3, CTaa = letters[1:3])))
  
  list_of_items.2 <- list(df1 = data.frame(x = 1:3, CTaa = letters[1:3]), 
                        df2 = data.frame(x = 6:9, CTaa = letters[8:11]))
  expect_equal(.checkBlanks(list_of_items.2, "CTaa"), 
               list(df1 = data.frame(x = 1:3, CTaa = letters[1:3]), 
                    df2 = data.frame(x = 6:9, CTaa = letters[8:11])))
  
})


test_that(".short.check works", {
  list_of_items <- list(df1 = data.frame(x = 1:3, CTaa = letters[1:3]), 
                        df2 = data.frame(x = 6:9, CTaa = letters[8:11]))
  expect_equal(.short.check(list_of_items, "CTaa"), 3)
  
})


test_that(".get.coord works", {
  data("scRep_example")
  expect_equal(.get.coord(scRep_example, "umap"),
               getdata("utils", "get.coord_data"))
})


# Test cases
test_that(".modifyBarcodes without ID works correctly", {
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
test_that(".modifyBarcodes with ID works correctly", {
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

test_that(".assignCT without ID works correctly", {
    test_contigs <- getCombined()[[1]][,c("barcode", "TCR1", "TCR2", "cdr3_nt1", "cdr3_nt2", "cdr3_aa1", "cdr3_aa2")]
    expect_equal(.assignCT("T",test_contigs),
                 getdata("utils", "assignCT_Tcell_data"))
})
#TODO Add B cell data to assignCT check
 

test_that(".colorizer works", {
  expect_identical(.colorizer(palette = "inferno", n= 5), 
                   c("#040404", "#611163", "#C53270", "#F69422", "#FFFE9E"))
  expect_identical(rev(.colorizer(palette = "Greens", n= 5)), 
                   c("#F6FBF4", "#CAE8C1", "#81C07A", "#30893B", "#004616"))
  
})



test_that(".theCall works", {
  expect_equal(.theCall(NULL, "aa", check.df = FALSE), "CTaa")
  expect_equal(.theCall(NULL, "nt", check.df = FALSE), "CTnt")
  expect_equal(.theCall(NULL, "genes", check.df = FALSE), "CTgene")
  expect_equal(.theCall(NULL, "strict", check.df = FALSE), "CTstrict")
})
#TODO .theCall Add custom header

test_that(".constructConDfAndParseTCR works", {
  expect_identical(
    .constructConDfAndParseTCR(
      getdata("utils", "constructConDfAndParseTCRInput")
    ),
    getdata("utils", "expected_con_df")
  )

# TODO: add more cases! This is not comprehensive.
# there was a case in the past where this passed but
# the function caused a segmentation fault.
})



test_that("makeGenes works for cellType T", {
    expect_identical(
        .makeGenes("T", getdata("utils", "makeGenes_T_input")),
        getdata("utils", "makeGenes_T_expected")
    )
})
# TODO makesGenes (cellType B)



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
