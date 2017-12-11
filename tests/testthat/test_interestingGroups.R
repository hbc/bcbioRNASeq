context("interestingGroups")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("accessor", {
    expect_equal(
        interestingGroups(bcb),
        "group"
    )
})

test_that("assignment", {
    expect_silent(
        interestingGroups(bcb) <- "sampleName"
    )
    # Interesting group must be in the metadata
    expect_error(
        interestingGroups(bcb) <- "XXX",
        "Interesting groups not defined in metadata: XXX"
    )

})
