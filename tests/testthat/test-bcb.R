context("bcb example dataset")

## Check to make sure tximport loaded correctly.
test_that("Assays", {
    assaySum <- function(x) {
        ## Handling NA values in aligned counts here.
        x[is.na(x)] <- 0L
        as.integer(round(sum(x)))
    }
    ## Transcript-to-gene counts aren't integer, so we're rounding here to
    ## check the values more easily.
    expect_identical(
        object = lapply(assays(bcb), assaySum),
        expected = list(
            counts = 179646L,
            tpm = 3456L,
            avgTxLength = 1389741,
            aligned = 138365L,
            normalized = 163383L,
            vst = 4044L,
            fpkm = 511409L
        )
    )
})
