test_that("bcbioRNASeq to DESeqDataSet", {
    x <- as(object, "DESeqDataSet")
    expect_s4_class(x, "DESeqDataSet")
    expect_identical(assayNames(x), "counts")
    expect_true(is.integer(assay(x)))
    expect_identical(
        object = colSums(assay(x)),
        expected = c(
            ## nolint start
            control_rep1 = 17661,
            control_rep2 = 60247,
            control_rep3 = 16106,
            fa_day7_rep1 = 29483,
            fa_day7_rep2 = 26270,
            fa_day7_rep3 = 29885
            ## nolint end
        )
    )
})

test_that("bcbioRNASeq to DESeqTransform", {
    x <- as(object, "DESeqTransform")
    expect_s4_class(x, "DESeqTransform")
})

test_that("bcbioRNASeq to DGEList (edgeR)", {
    x <- as(object, "DGEList")
    expect_s4_class(x, "DGEList")
    expect_true(validObject(x))
    counts <- x[["counts"]]
    expect_is(counts, "matrix")
    expect_false(is.integer(counts))
    expect_identical(
        object = round(colSums(counts)),
        expected = c(
            ## nolint start
            control_rep1 = 17659,
            control_rep2 = 60245,
            control_rep3 = 16105,
            fa_day7_rep1 = 29482,
            fa_day7_rep2 = 26272,
            fa_day7_rep3 = 29883
            ## nolint end
        )
    )
})
