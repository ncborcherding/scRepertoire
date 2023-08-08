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

# off.the.chain
# checkBlanks
# groupList

test_that("checkList works", {
    data("contig_list")
    expect_identical(checkList(contig_list), contig_list)
    expect_identical(checkList(contig_list[[1]]), contig_list[[1]])
    # no idea what to put to make the stop message happen. perhaps with data.table(contig_list) ? but then it shouldnt fail :P
})

test_that("checkContigs works", {
    input <- list(
        df1 = data.frame(a = c("x", "", "z"), b = c("1", "2", "3")),
        df2 = data.frame(c = c("foo", "bar", ""), d = c("", "spam", "eggs"))
    )
    expected <- checkContigs(input)
    expect_true(is.list(expected))
    expect_true(is.data.frame(expected[[1]]))
    expect_equal(expected[[1]]$a, c("x", NA, "z"))
    expect_equal(expected[[1]]$b, c("1", "2", "3"))
    expect_true(is.data.frame(expected[[2]]))
    expect_equal(expected[[2]]$c, c("foo", "bar", NA))
    expect_equal(expected[[2]]$d, c(NA, "spam", "eggs"))
})

#bound.input.return
#get.coord
#checkSingleObject
#grabMeta
#modifyBarcodes
#removingNA
#removingMulti

test_that("filteringMulti works", {
    source("testdata/utils-testdata.R")
    expect_identical(filteringMulti(head(contig_list[[1]])), filteringMulti_expected)
})

# TODO: tests for the rest of the functions in R/utils.R
