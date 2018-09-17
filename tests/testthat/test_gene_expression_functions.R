context("Gene Expression Functions")

gene2symbol <- gene2symbol(bcb_small)
geneIDs <- head(gene2symbol[["geneID"]])
geneNames <- head(gene2symbol[["geneName"]])



# plotGenderMarkers ============================================================
with_parameters_test_that(
    "plotGenderMarkers", {
        x <- plotGenderMarkers(object)
        expect_is(x, "ggplot")
    },
    object = list(
        bcbioRNASeq = bcb_small,
        DESeqDataSet = dds_small,
        DESeqTransform = vst_small
    )
)



# plotGene =====================================================================
# FIXME parameterize
test_that("plotGene : bcbioRNASeq", {
    # facet
    p <- plotGene(
        object = bcb_small,
        genes = geneNames,
        normalized = "vst",
        interestingGroups = "sampleName",
        style = "facet"
    )
    expect_is(p, "ggplot")

    # wide
    p <- plotGene(
        object = bcb_small,
        genes = geneNames,
        normalized = "tpm",
        style = "wide"
    )
    expect_is(p, "ggplot")
})

test_that("plotGene : DESeqDataSet", {
    p <- plotGene(dds_small, genes = geneNames)
    expect_is(p, "ggplot")
})

test_that("plotGene : DESeqTransform", {
    # rlog
    p <- plotGene(rlog_small, genes = geneNames)
    expect_is(p, "ggplot")

    # vst
    p <- plotGene(vst_small, genes = geneNames)
    expect_is(p, "ggplot")
})



# plotHeatmap ==================================================================
# FIXME parameterize
test_that("plotHeatmap : bcbioRNASeq", {
    genes <- head(rownames(bcb_small), n = 100L)
    p <- plotHeatmap(bcb_small[genes, ])
    expect_is(p, "pheatmap")
})

test_that("plotHeatmap : DESeqDataSet", {
    genes <- head(rownames(dds_small), n = 20L)
    p <- plotHeatmap(dds_small[genes, ])
    expect_is(p, "pheatmap")
})
