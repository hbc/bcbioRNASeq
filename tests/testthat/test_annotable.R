context("annotable")

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
        c(ensgene = "character",
          symbol = "character",
          description = "character",
          biotype = "character",
          broadClass = "character")
    )
})
