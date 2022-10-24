test_that("bcbioRNASeq", {
    ## Fast mode skips import of aligned counts.
    bcb <- bcbioRNASeq(uploadDir, fast = TRUE)
    expect_false("aligned" %in% assayNames(bcb))
    bcb <- slotAlignedCounts(bcb)
    expect_true("aligned" %in% assayNames(bcb))
    expect_error(
        object = slotAlignedCounts(bcb),
        regexp = "aligned"
    )
})
