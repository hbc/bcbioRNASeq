context("data : bcb")

with_parameters_test_that(
    "Assays", {
        expect_is(object, "matrix")
        # Check to make sure tximport loaded correctly.
        # Transcript-to-gene counts aren't integer, so we're rounding here to
        # check the values more easily.
        expect_identical(
            object = round(sum(object)),
            expected = sum
        )
    },
    object = as.list(assays(bcb)),
    # nolint start
    sum = list(
        counts = 143758911,
        tpm = 5249973,
        avgTxLength = 2696963,
        normalized = 137992370,
        vst = 32371,
        fpkm = 3800855
    )
    # nolint end
)
