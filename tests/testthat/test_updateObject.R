context("updateObject")

test_that("v0.1.4", {
    loadRemoteData("http://bcbiornaseq.seq.cloud/v0.1.4/bcb.rda", quiet = TRUE)
    expect_identical(
        metadata(bcb)[["version"]],
        package_version("0.1.4")
    )
    updated <- suppressMessages(updateObject(bcb))
    expect_identical(
        metadata(updated)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(updated)[["previousVersion"]],
        package_version("0.1.4")
    )
})
