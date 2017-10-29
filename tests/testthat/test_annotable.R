context("annotable")

test_that("annotable", {
    anno <- annotable(bcb)
    expect_true(is.data.frame(anno))
    expect_equal(
        colnames(anno),
        c("ensgene",
          "symbol",
          "description",
          "biotype",
          "broadClass")
    )
    expect_equal(
        rownames(anno)[1:3],
        c("ENSMUSG00000000001",
          "ENSMUSG00000000003",
          "ENSMUSG00000000028")
    )
})
