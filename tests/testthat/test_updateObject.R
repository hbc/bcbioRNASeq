context("updateObject")

test_that("v0.1.4", {
    expect_error(validObject(bcb_invalid))
    expect_identical(
        metadata(bcb_invalid)[["version"]],
        package_version("0.1.4")
    )
    organism <- metadata(bcb_invalid)[["organism"]]
    rowRanges <- genes(organism)
    bcb <- suppressWarnings(updateObject(bcb_invalid, rowRanges = rowRanges))
    expect_identical(
        metadata(bcb)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(bcb)[["previousVersion"]],
        package_version("0.1.4")
    )
})
