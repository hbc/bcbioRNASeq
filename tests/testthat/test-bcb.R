## Check to make sure tximport loaded correctly.
test_that("Assays", {
    assaySum <- function(x) {
        ## Handling NA values in aligned counts here.
        x[is.na(x)] <- 0L
        as.integer(round(sum(x)))
    }
    ## Transcript-to-gene counts aren't integer, so we're rounding here to
    ## check the values more easily.
    sums <- unlist(lapply(assays(bcb), assaySum))
    sums <- sums[sort(names(sums))]
    expect_identical(
        object = sums,
        expected = c(
            aligned = 138365L,
            avgTxLength = 1389741L,
            counts = 179646L,
            fpkm = 511409L,
            normalized = 163383L,
            tpm = 3456L,
            vst = 4044L
        )
    )
})
