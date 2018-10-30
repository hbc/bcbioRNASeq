context("Differential Expression")

data(bcb, deseq, envir = environment())
dds_small <- as(deseq, "DESeqDataSet")
vst_small <- as(deseq, "DESeqTransform")
res_small <- as(deseq, "DESeqResults")

g2s <- Gene2Symbol(bcb)
geneIDs <- head(g2s[["geneID"]])
geneNames <- head(g2s[["geneName"]])



# alphaSummary =================================================================
test_that("alphaSummary : DESeqDataSet", {
    object <- dds_small

    # Default, no contrast specified.
    x <- alphaSummary(object)
    expect_is(x, "matrix")
    expect_equal(
        object = x,
        expected = matrix(
            # nolint start
            data = c(
                115,  88,  58,  31,  19,
                139, 120, 100,  59,  14,
                  6,   6,   6,    6,  6,
                  0,   0,   0,    0,  39
            ),
            # nolint end
            nrow = 4L,
            ncol = 5L,
            byrow = TRUE,
            dimnames = list(
                c(
                    "LFC > 0 (up)",
                    "LFC < 0 (down)",
                    "outliers [1]",
                    "low counts [2]"
                ),
                c(1e-01, 0.5e-01, 1e-02, 1e-03, 1e-06)
            )
        )
    )

    # Contrast vector.
    expect_identical(
        object = alphaSummary(
            object = object,
            contrast = c("treatment", "folic_acid", "control")
        ),
        expected = x
    )

    # Contrast name.
    expect_identical(
        object = alphaSummary(
            object = object,
            name = "treatment_folic_acid_vs_control"
        ),
        expected = x
    )
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
        DESeqAnalysis = list(object = deseq)
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
        DESeqAnalysis = list(object = deseq)
    )
)



# plotMA =======================================================================
with_parameters_test_that(
    "plotMA", {
        x <- plotMA(object)
        expect_is(x, "ggplot")

        # Check geom classes.
        geomtype <- vapply(
            X = x[["layers"]],
            FUN = function(x) {
                class(x[["geom"]])[[1L]]
            },
            FUN.VALUE = character(1L)
        )
        expect_identical(
            object = geomtype,
            expected = c("GeomHline", "GeomPoint", "GeomLogticks")
        )

        # Check plot labels.
        expect_identical(
            object = x[["labels"]][["y"]],
            expected = "log2 fold change"
        )
        expect_identical(
            object = x[["labels"]][["x"]],
            expected = "mean expression across all samples"
        )

        # Directional support.
        x <- plotMA(
            object = object,
            direction = "up",
            sigPointColor = "red"
        )
        expect_is(x, "ggplot")
        x <- plotMA(
            object = object,
            direction = "down",
            sigPointColor = "green"
        )
        expect_is(x, "ggplot")

        # Label the top genes.
        args <- list(object = object, ntop = 10L)
        if (!is(object, "DESeqAnalysis")) {
            args[["gene2symbol"]] <- g2s
        }
        x <- do.call(what = plotMA, args = args)
        expect_is(x, "ggplot")

        # Label specific genes.
        args <- list(object = object, genes = geneNames)
        if (!is(object, "DESeqAnalysis")) {
            args[["gene2symbol"]] <- g2s
        }
        x <- do.call(what = plotMA, args = args)
        expect_is(x, "ggplot")

        # DataFrame return.
        x <- plotMA(object, return = "DataFrame")
        expect_s4_class(x, "DataFrame")
        expect_true("isDE" %in% colnames(x))
    },
    object = list(
        DESeqAnalysis = deseq,
        DESeqResults = res_small
    )
)



# plotVolcano ==================================================================
with_parameters_test_that(
    "plotVolcano", {
        x <- plotVolcano(object)
        expect_is(x, "ggplot")

        # Enable histograms.
        x <- plotVolcano(object, histograms = TRUE)
        expect_is(x, "ggplot")

        # Directional support.
        x <- plotVolcano(
            object = object,
            direction = "up",
            sigPointColor = "red"
        )
        expect_is(x, "ggplot")
        x <- plotVolcano(
            object = object,
            direction = "down",
            sigPointColor = "green"
        )
        expect_is(x, "ggplot")

        # Label the top genes.
        args <- list(object = object, ntop = 5L)
        if (!is(object, "DESeqAnalysis")) {
            args[["gene2symbol"]] <- g2s
        }
        x <- do.call(what = plotVolcano, args = args)
        expect_is(x, "ggplot")

        # Label specific genes.
        args <- list(object = object, genes = geneNames)
        if (!is(object, "DESeqAnalysis")) {
            args[["gene2symbol"]] <- g2s
        }
        x <- do.call(what = plotVolcano, args = args)
        expect_is(x, "ggplot")

        # Return DataFrame.
        x <- plotVolcano(object, return = "DataFrame")
        expect_s4_class(x, "DataFrame")
    },
    object = list(
        DESeqAnalysis = deseq,
        DESeqResults = res_small
    )
)



# FIXME Move to generators test file...
# DESeqResultsTables ===========================================================
with_parameters_test_that(
    "DESeqResultsTables", {
        expect_s4_class(
            object = DESeqResultsTables(object),
            class = "DESeqResultsTables"
        )
    },
    object = list(
        DESeqAnalysis = deseq,
        DESeqResults = res_small
    )
)
