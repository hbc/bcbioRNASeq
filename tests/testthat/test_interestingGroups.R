context("interestingGroups")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("Accessor", {
    expect_equal(
        interestingGroups(bcb),
        "group"
    )
})

test_that("Assignment", {
    expect_silent(
        interestingGroups(bcb) <- "sampleName"
    )
    # Interesting group must be in the metadata
    expect_error(
        interestingGroups(bcb) <- "XXX",
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(x\\)."
        )
    )
})
