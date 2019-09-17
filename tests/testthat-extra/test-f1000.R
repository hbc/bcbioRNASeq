## nolint start



object <- import("https://github.com/hbc/bcbioRNASeq/raw/f1000v2/data/bcb.rda")
object <- updateObject(object)
expect_s4_class(object, "bcbioRNASeq")



context("F1000 workflow paper")

test_that("counts", {
    raw <- counts(object, normalized = FALSE)
    expect_is(raw, "matrix")

    normalized <- counts(object, normalized = TRUE)
    expect_is(normalized, "matrix")

    tpm <- counts(object, normalized = "tpm")
    expect_is(tpm, "matrix")

    rlog <- counts(object, normalized = "rlog")
    expect_is(rlog, "matrix")

    vst <- counts(object, normalized = "vst")
    expect_is(vst, "matrix")
})

test_that("Quality control", {
    p <- plotTotalReads(object)
    expect_s3_class(p, "ggplot")

    p <- plotMappingRate(object)
    expect_s3_class(p, "ggplot")

    p <- plotExonicMappingRate(object)
    expect_s3_class(p, "ggplot")

    p <- plotIntronicMappingRate(object)
    expect_s3_class(p, "ggplot")

    expect_warning(
        object = plotGenesDetected(object),
        regexp = "deprecated"
    )

    p <- plotGeneSaturation(object)
    expect_s3_class(p, "ggplot")

    expect_warning(
        object = plotCountsPerGene(object),
        regexp = "deprecated"
    )

    p <- plotCountsPerFeature(object, geom = "density")
    expect_identical(nrow(p$data), 457704L)
    expect_match(p$labels$subtitle, "38142")

    p <- plotCountDensity(object)
    expect_s3_class(p, "ggplot")

    ## > p <- plotMeanSD(object)
    ## > expect_s3_class(p, "ggplot")

    ## > p <- plotDispEsts(object)
    ## > expect_is(p, "list")

    p <- plotCorrelationHeatmap(object)
    expect_s3_class(p, "pheatmap")

    p <- plotPCA(object)
    expect_s3_class(p, "ggplot")

    p <- plotPCACovariates(object, fdr = 0.1)
    expect_is(p, "list")
})

test_that("Differential expression", {
    results <- DESeq2::results

    dds <- as(object, "DESeqDataSet")
    design(dds) <- ~ day
    dds <- DESeq(dds)
    vst <- varianceStabilizingTransformation(dds)

    x <- capture.output(
        alphaSummary(
            object = dds,
            contrast = c(
                factor = "day",
                numerator = "7",
                denominator = "0"
            ),
            alpha = c(0.1, 0.05)
        )
    )
    expect_identical(
        object = x[[2L]],
        expected = "LFC > 0 (up)    1521  1102"
    )

    res <- results(object = dds, name = "day_7_vs_0", alpha = 0.05)
    expect_s4_class(res, "DESeqResults")

    p <- plotMeanAverage(res)
    expect_s3_class(p, "ggplot")

    p <- plotVolcano(res)
    expect_s3_class(p, "ggplot")

    ## Note that the generic formals have been renamed here.
    ## object : results
    ## DESeqTransform : counts
    p <- plotDEGHeatmap(
        object = res,
        DESeqTransform = vst
    )
    expect_s3_class(p, "pheatmap")

    ## lfc argument renamed to lfcThreshold.
    res_tbl <- resultsTables(res, lfcThreshold = 1)
    expect_is(res_tbl, "list")
    expect_is(res_tbl[[1L]], "tbl_df")

    expect_output(topTables(res_tbl, n = 5))
})



context("R Markdown templates")

## FIXME



unlink(dir, recursive = TRUE)



## nolint end
