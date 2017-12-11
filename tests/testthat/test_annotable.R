context("annotable")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))

test_that("annotable", {
    anno <- annotable(bcb)
    expect_is(anno, "data.frame")
    expect_equal(
        nrow(anno),
        nrow(bcb)
    )
    expect_equal(
        rownames(anno),
        rownames(bcb)
    )
    expect_equal(
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
