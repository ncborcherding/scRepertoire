test_that("'%!in%' works", {
    vector <- c(1, 2, 3, 4, 5)
    
    expect_true(0 %!in% vector)
    expect_true(6 %!in% vector)
    expect_false(3 %!in% vector)
    expect_false(5 %!in% vector)
    expect_true(1 %!in% NULL)
    expect_true(list(1) %!in% NA)
})
