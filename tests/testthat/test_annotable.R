context("annotable")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("annotable", {
    anno <- annotable(bcb)
    expect_is(anno, "data.frame")
    expect_identical(
        nrow(anno),
        nrow(bcb)
    )
    expect_identical(
        rownames(anno),
        rownames(bcb)
    )
    expect_identical(
        lapply(anno, class),
        list(ensgene = "character",
             symbol = "character",
             description = "character",
             biotype = "character",
             broadClass = "character",
             geneSeqStart = "integer",
             geneSeqEnd = "integer",
             seqName = "character",
             seqStrand = "integer",
             seqCoordSystem = "character",
             # This is returning `AsIs` instead of `list`
             entrez = "AsIs")
    )
})
