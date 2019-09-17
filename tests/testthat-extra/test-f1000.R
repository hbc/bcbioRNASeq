## nolint start



context("F1000 workflow paper")

library(DESeq2)
library(DEGreport)
library(DESeqAnalysis)

dir <- "f1000"
unlink(dir, recursive = TRUE)

loadRemoteData("https://github.com/hbc/bcbioRNASeq/raw/f1000v2/data/bcb.rda")
bcb <- updateObject(bcb)
expect_s4_class(bcb, "bcbioRNASeq")

raw <- counts(bcb, normalized = FALSE)
expect_is(raw, "matrix")

normalized <- counts(bcb, normalized = TRUE)
expect_is(normalized, "matrix")

tpm <- counts(bcb, normalized = "tpm")
expect_is(tpm, "matrix")

rlog <- counts(bcb, normalized = "rlog")
expect_is(rlog, "matrix")

vst <- counts(bcb, normalized = "vst")
expect_is(vst, "matrix")

p <- plotTotalReads(bcb)
expect_s3_class(p, "ggplot")

p <- plotMappingRate(bcb)
expect_s3_class(p, "ggplot")

p <- plotExonicMappingRate(bcb)
expect_s3_class(p, "ggplot")

p <- plotIntronicMappingRate(bcb)
expect_s3_class(p, "ggplot")

expect_warning(
    object = plotGenesDetected(bcb),
    regexp = "deprecated"
)

p <- plotGeneSaturation(bcb)
expect_s3_class(p, "ggplot")

expect_warning(
    object = plotCountsPerGene(bcb),
    regexp = "deprecated"
)

p <- plotCountsPerFeature(bcb, geom = "density")
expect_identical(nrow(p$data), 457704L)
expect_match(p$labels$subtitle, "38142")

p <- plotCountDensity(bcb)
expect_s3_class(p, "ggplot")

p <- plotMeanSD(bcb)
expect_s3_class(p, "ggplot")

p <- plotDispEsts(bcb)
expect_is(p, "list")

p <- plotCorrelationHeatmap(bcb)
expect_s3_class(p, "pheatmap")

p <- plotPCA(bcb)
expect_s3_class(p, "ggplot")

p <- plotPCACovariates(bcb, fdr = 0.1)
expect_is(p, "list")

dds <- as(bcb, "DESeqDataSet")
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

## FIXME This doesn't work with SimpleDataFrameList (see above).
## FIXME Need to add `list` method back to support F1000 paper.
topTables(res_tbl, n = 5)

unlink(dir, recursive = TRUE)



## nolint end
