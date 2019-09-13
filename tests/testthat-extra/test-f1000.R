## nolint start



context("F1000 workflow paper")

## FIXME Switch paper to use `deg()` instead of `significants()`.

require(DEGreport)
require(DESeqAnalysis)
require(clusterProfiler)

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

p <- plotGenesDetected(bcb)
expect_s3_class(p, "ggplot")

p <- plotGeneSaturation(bcb)
expect_s3_class(p, "ggplot")

p <- plotCountsPerGene(bcb)
expect_s3_class(p, "ggplot")

## FIXME There's a difference in n now: 38142 vs. 27041.
## > p <- plotCountsPerFeature(bcb, geom = "density")

## FIXME This is in paper, but now defunct in package.
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
    x[[2L]],
    "LFC > 0 (up)    1521  1102"
)

res <- results(
    object = dds,
    name = "day_7_vs_0",
    alpha = 0.05
)
expect_s4_class(res, "DESeqResults")

p <- plotMeanAverage(res)
expect_s3_class(p, "ggplot")

p <- plotVolcano(res)
expect_s3_class(p, "ggplot")

## Note that the generic formals have been renamed here.
p <- plotDEGHeatmap(
    object = res,
    DESeqTransform = vst
)
expect_s3_class(p, "pheatmap")

## Note that "log" argument has been changed to "log2" and flipped.
p <- degPlot(
    bcb,
    xs = "day",
    res = res,
    n = 3,
    slot = "vst",
    log2 = FALSE,
    ann = c("geneID", "geneName"),
)
expect_s3_class(p, "ggplot")

dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
res_lrt <- results(dds_lrt)

ma <- counts(bcb, "vst")[significants(res, fc = 2), ]
expect_is(ma, "matrix")
expect_identical(dim(ma), c(1036L, 12L))

res_patterns <- degPatterns(
    ma = ma,
    metadata = colData(bcb),
    time = "day",
    minc = 60
)
expect_is(res_patterns, "list")
## FIXME Take this code out in paper...
## > res_patterns[["plot"]]
expect_s3_class(res_patterns[["plot"]], "ggplot")

x <- subset(res_patterns[["df"]], cluster == 3, select = "genes")
expect_is(x, "data.frame")

## FIXME Return has changed from tibble list to SplitDataFrameList.
## FIXME Change the default back in DESeqAnalysis?
res_tbl <- resultsTables(res, lfc = 1)
expect_is(res_tbl, "list")
expect_is(res_tbl[[1L]], "tbl_df")

## FIXME This doesn't work with SimpleDataFrameList (see above).
## FIXME Need to add `list` method back to support F1000 paper.
topTables(res_tbl, n = 5)

## Generate the gene IDs for the DE list of genes and the background list of
## genes.
all_genes <- rownames(res) %>%
    .[!is.na(res[["padj"]])] %>%
    as.character()
sig_genes <- significants(
    object = res,
    fc = 1,
    padj = 0.05
)

## Generate fold change values for significant results.
sig_results <- as.data.frame(res)[sig_genes, ]
fold_changes <- sig_results$log2FoldChange
names(fold_changes) <- rownames(sig_results)

## Run GO enrichment analysis.

ego <- enrichGO(
    gene = sig_genes,
    OrgDb = "org.Mm.eg.db",
    keyType = "ENSEMBL",
    ont = "BP",
    universe = all_genes,
    qvalueCutoff = 0.05,
    readable = TRUE
)

ego_summary <- slot(ego, "result") %>%
    as_tibble() %>%
    camel()

# Dotplot of top 25

## FIXME wrong orderBy parameter; set to default `orderBy = "x"`
## Report to DOSE package author.
dotplot(ego, showCategory = 25)

## Enrichment plot of top 25.
## FIXME This is no longer in clusterProfiler...
enrichMap(ego, n = 25, vertex.label.cex = 0.5)



## nolint end
