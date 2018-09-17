context("Differential Expression Functions")

gene2symbol <- gene2symbol(bcb_small)
geneIDs <- head(gene2symbol[["geneID"]])
geneNames <- head(gene2symbol[["geneName"]])



# alphaSummary =================================================================
test_that("alphaSummary : DESeqDataSet", {
    # Default, no contrast specified.
    object <- alphaSummary(dds_small)
    expect_is(object, "knitr_kable")
    expect_true(grepl("1e-06", object[[1L]]))

    # Contrast vector.
    object <- alphaSummary(
        object = dds_small,
        contrast = c("treatment", "folic_acid", "control")
    )
    expect_is(object, "knitr_kable")

    # Contrast name.
    object <- alphaSummary(
        object = dds_small,
        name = "treatment_folic_acid_vs_control"
    )
    expect_is(object, "knitr_kable")
})



# plotDEGHeatmap ===============================================================
with_parameters_test_that(
    "plotDEGHeatmap", {
        expect_is(
            object = do.call(what = plotDEGHeatmap, args = args),
            class = "pheatmap"
        )
    },
    args = list(
        DESeqResults = list(
            object = res_small,
            counts = vst_small
        ),
        DESeqAnalysis = list(object = deseq_small)
    )
)



# plotDEGPCA ===================================================================
with_parameters_test_that(
    "plotDEGPCA", {
        expect_is(
            object = do.call(what = plotDEGPCA, args = args),
            class = "ggplot"
        )
    },
    args = list(
        DESeqResults = list(
            object = res_small,
            counts = vst_small
        ),
        DESeqAnalysis = list(object = deseq_small)
    )
)



# plotMA =======================================================================
test_that("plotMA : DESeqResults", {
    object <- plotMA(res_small)
    expect_is(object, "ggplot")

    # Check geom classes
    geomtype <- vapply(
        X = object[["layers"]],
        FUN = function(object) {
            class(object[["geom"]])[[1L]]
        },
        FUN.VALUE = character(1L)
    )
    expect_identical(
        geomtype,
        c("GeomHline", "GeomPoint", "GeomLogticks")
    )

    # Check plot labels
    expect_identical(
        object[["labels"]][["y"]],
        "log2 fold change"
    )
    expect_identical(
        object[["labels"]][["x"]],
        "mean expression across all samples"
    )
})

test_that("plotMA : Specific genes", {
    object <- plotMA(
        object = res_small,
        genes = geneNames,
        gene2symbol = gene2symbol
    )
    expect_is(object, "ggplot")
})

test_that("plotMA: ntop mode", {
    object <- plotMA(
        object = res_small,
        ntop = 10L,
        gene2symbol = gene2symbol
    )
    expect_is(object, "ggplot")
})

test_that("plotMA : Directional support", {
    # Upregulated
    object <- plotMA(
        object = res_small,
        direction = "up",
        sigPointColor = "red"
    )
    expect_is(object, "ggplot")

    # Downregulated
    object <- plotMA(
        object = res_small,
        direction = "down",
        sigPointColor = "green"
    )
    expect_is(object, "ggplot")
})

test_that("plotMA : DataFrame return", {
    object <- plotMA(res_small, return = "DataFrame")
    expect_s4_class(object, "DataFrame")
    expect_true("isDE" %in% colnames(object))
})



# plotVolcano ==================================================================
test_that("plotVolcano : DESeqResults", {
    object <- plotVolcano(res_small, gene2symbol = gene2symbol)
    expect_is(object, "ggplot")

    # Enable histograms
    object <- plotVolcano(res_small, histograms = TRUE)
    expect_is(object, "ggplot")

    # Label the top genes
    object <- plotVolcano(res_small, ntop = 5L, gene2symbol = gene2symbol)
    expect_is(object, "ggplot")

    # Label specific genes
    object <- plotVolcano(
        object = res_small,
        genes = geneNames,
        gene2symbol = gene2symbol
    )
    expect_is(object, "ggplot")

    # Directional support
    object <- plotVolcano(
        object = res_small,
        direction = "up",
        sigPointColor = "red"
    )
    expect_is(object, "ggplot")
    object <- plotVolcano(
        object = res_small,
        direction = "down",
        sigPointColor = "green"
    )
    expect_is(object, "ggplot")

    # Return DataFrame
    object <- plotVolcano(res_small, return = "DataFrame")
    expect_s4_class(object, "DataFrame")
})



# resultsTables ================================================================
with_parameters_test_that(
    "resultsTables", {
        expect_s4_class(
            object = resultsTables(object),
            class = "DESeqResultsTables"
        )
    },
    object = list(
        DESeqAnalysis = deseq_small,
        DESeqResults = res_small
    )
)
