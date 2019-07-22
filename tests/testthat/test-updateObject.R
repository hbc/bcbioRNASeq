context("updateObject")

## Load a legacy object that doesn't contain rowRanges.
load("bcb_invalid.rda")
object <- bcb_invalid

test_that("updateObject", {
    expect_error(
        object = validObject(object),
        regexp = "rowRanges"
    )
    expect_identical(
        object = slot(object, "metadata")[["version"]],
        expected = package_version("0.1.4")
    )

    x <- suppressWarnings(updateObject(object))
    expect_s4_class(x, "bcbioRNASeq")
    expect_true(validObject(x))
    expect_identical(
        object = metadata(x)[["version"]],
        expected = packageVersion
    )
    expect_identical(
        object = metadata(x)[["previousVersion"]],
        expected = package_version("0.1.4")
    )
})
