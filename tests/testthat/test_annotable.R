context("annotable")

test_that("annotable", {
    x <- annotable(bcb)
    expect_is(x, "data.frame")
    expect_identical(
        nrow(x),
        nrow(x)
    )
    expect_identical(
        rownames(x),
        rownames(bcb)
    )
    expect_identical(
        lapply(x, class),
        list(
            ensgene = "character",
            symbol = "character",
            description = "character",
            biotype = "character",
            broadClass = "character",
            geneSeqStart = "integer",
            geneSeqEnd = "integer",
            seqName = "character",
            seqStrand = "integer",
            seqCoordSystem = "character",
            # FIXME Ensure `list` return instead of `AsIs`
            entrez = "AsIs")
    )
})
