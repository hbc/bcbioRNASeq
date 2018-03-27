context("updateObject")

test_that("updateObject", {
    invalid <- bcb_small
    metadata(invalid)[["version"]] <- package_version("0.1.7")
    expect_error(
        validObject(invalid),
        "metadata\\(object\\)\\[\\[\"version\"\\]\\] >= 0.2 is not TRUE"
    )
    organism <- metadata(invalid)[["organism"]]
    rowRanges <- genes(organism)
    valid <- suppressWarnings(updateObject(invalid, rowRanges = rowRanges))
    expect_identical(
        metadata(valid)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(valid)[["previousVersion"]],
        package_version("0.1.7")
    )
})
