quietedTCRvec <- getdata("quietVDJ", "quietedTCRvector")
quietedBCRvec <- getdata("quietVDJ", "quietedBCRvector")

test_that("quietTCRgenes works", {
    expect_setequal(quietTCRgenes(quietedTCRvec), quietedTCRvec)
})

test_that("quietBCRgenes works", {
    expect_setequal(quietBCRgenes(quietedBCRvec), quietedBCRvec)
})

test_that("quietVDJgenes works", {
    expect_setequal(
        quietVDJgenes(c(quietedTCRvec, quietedBCRvec)),
        quietVDJgenes(union(quietedTCRvec, quietedBCRvec))
    )
    expect_setequal(
        quietVDJgenes(setdiff(quietedTCRvec, quietedBCRvec)),
        character(0)
    )
})
