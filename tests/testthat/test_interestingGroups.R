context("interestingGroups")

test_that("Accessor", {
    expect_identical(
        interestingGroups(bcb_small),
        "treatment"
    )
})

test_that("Assignment", {
    expect_silent(
        interestingGroups(bcb_small) <- "sampleName"
    )
    # Interesting group must be in the metadata
    expect_error(
        interestingGroups(bcb_small) <- "XXX",
        paste(
            "is_subset :",
            "The element 'XXX' in interestingGroups is not in",
            "colnames\\(x\\)."
        )
    )
})
