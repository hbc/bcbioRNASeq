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

    # contrast vector
    object <- alphaSummary(
        object = dds_small,
        contrast = c("treatment", "folic_acid", "control")
    )
    expect_is(object, "knitr_kable")

    # contrast name
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
            object = plotDEGHeatmap(
                results = res_small,
                counts = counts
            ),
            class = "pheatmap"
        )
    },
    counts = list(
        bcbioRNASeq = bcb_small,
        DESeqDataSet = dds_small,
        DESeqTransform = vst_small
    )
)

test_that("plotDEGHeatmap : No DEGs", {
    args <- list(
        results = res_small,
        counts = bcb_small,
        lfcThreshold = Inf
    )
    expect_null(do.call(what = plotDEGHeatmap, args = args))
    expect_warning(
        object = do.call(what = plotDEGHeatmap, args = args),
        regexp = "No significant DEGs to plot"
    )
})



# plotDEGPCA ===================================================================
with_parameters_test_that(
    "plotDEGPCA", {
        expect_is(
            object = plotDEGPCA(
                results = res_small,
                counts = counts
            ),
            class = "ggplot"
        )
    },
    counts = list(
        bcbioRNASeq = bcb_small,
        DESeqDataSet = dds_small,
        DESeqTransform = vst_small
    )
)

test_that("plotDEGPCA : No DEGs", {
    args <- list(
        results = res_small,
        counts = bcb_small,
        lfcThreshold = Inf
    )
    expect_null(do.call(what = plotDEGPCA, args = args))
    expect_warning(
        object = do.call(what = plotDEGPCA, args = args),
        regexp = "No significant DEGs to plot"
    )
})



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
        object[["labels"]][["object"]],
        "mean expression across all samples"
    )
})

test_that("plotMA : Specific genes", {
    object <- plotMA(
        object = res_small,
        genes = genes,
        gene2symbol = gene2symbol(bcb_small)
    )
    expect_is(object, "ggplot")
})

test_that("plotMA: ntop mode", {
    object <- plotMA(
        object = res_small,
        ntop = 10L,
        gene2symbol = gene2symbol(bcb_small)
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
    object <- plotVolcano(res_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(object, "ggplot")

    # Directional support
    object <- plotVolcano(res_small, direction = "up", sigPointColor = "red")
    expect_is(object, "ggplot")
    object <- plotVolcano(res_small, direction = "down", sigPointColor = "green")
    expect_is(object, "ggplot")

    # Return DataFrame
    object <- plotVolcano(res_small, return = "DataFrame")
    expect_s4_class(object, "DataFrame")
})

# FIXME Improve counts pass-in checks
# FIXME Check gene symbols



# resultsTables ================================================================
test_that("resultsTables : Default return with local files only", {
    object <- resultsTables(
        results = res_small,
        counts = dds_small,
        lfcThreshold = lfc,
        summary = FALSE,
        write = FALSE
    )
    expect_identical(class(object), "list")
    expect_identical(
        lapply(object, class),
        list(
            deg = "data.frame",
            degLFC = "data.frame",
            degLFCUp = "data.frame",
            degLFCDown = "data.frame",
            all = "data.frame",
            contrast = "character",
            alpha = "numeric",
            lfcThreshold = "numeric"
        )
    )
    # Ensure that the counts columns are correct
    expect_identical(
        rownames(dds_small),
        rownames(object[["all"]])
    )
    expect_identical(
        counts(dds_small, normalized = TRUE),
        as.matrix(object[["all"]][, colnames(assay(dds_small))])
    )
})

test_that("resultsTables : Summary and write support", {
    # This is also capturing the contents of the list return. Not sure how to
    # fix here for knitr asis_output
    dir <- file.path(getwd(), "resultsTables")
    output <- capture.output(
        resultsTables(
            results = res_small,
            counts = dds_small,
            lfcThreshold = lfc,
            summary = TRUE,
            headerLevel = 2L,
            write = TRUE,
            dir = "resultsTables"
        ),
        type = "output"
    ) %>%
        # Just evaluate the R Markdown output
        head(28L)
    expect_identical(
        output,
        c(
            "",
            "",
            "## treatment folic acid vs control",
            "",
            "",
            "",
            "### Summary statistics",
            "",
            "",
            "- 500 genes in counts matrix",
            "- Base mean > 0: 500 genes (non-zero)",
            "- Base mean > 1: 500 genes",
            "- Alpha: 0.1",
            "- LFC threshold: 0.25",
            "- DEG pass alpha: 254 genes",
            "- DEG LFC up: 115 genes",
            "- DEG LFC down: 139 genes",
            "",
            "",
            "",
            "### Results tables",
            "",
            "",
            paste0(
                "- [`treatment_folic_acid_vs_control_all.csv.gz`]",
                "(",
                file.path(
                    dir,
                    "treatment_folic_acid_vs_control_all.csv.gz"
                ),
                "): ",
                "All genes, sorted by Ensembl identifier."
            ),
            paste0(
                "- [`treatment_folic_acid_vs_control_deg.csv.gz`]",
                "(",
                file.path(
                    dir,
                    "treatment_folic_acid_vs_control_deg.csv.gz"
                ),
                "): ",
                "Genes that pass the alpha (FDR) cutoff."
            ),
            paste0(
                "- [`treatment_folic_acid_vs_control_deg_lfc_up.csv.gz`]",
                "(",
                file.path(
                    dir,
                    "treatment_folic_acid_vs_control_deg_lfc_up.csv.gz"
                ),
                "): ",
                "Upregulated DEG; positive log2 fold change."
            ),
            paste0(
                "- [`treatment_folic_acid_vs_control_deg_lfc_down.csv.gz`]",
                "(",
                file.path(
                    dir,
                    "treatment_folic_acid_vs_control_deg_lfc_down.csv.gz"
                ),
                "): ",
                "Downregulated DEG; negative log2 fold change."
            ),
            ""
        )
    )
})

# Providing a corresponding DESeqDataSet for counts is recommended
test_that("resultsTables : DESeqResults minimal mode", {
    object <- resultsTables(results = res_small)
    expect_is(object, "list")
})

if (file.exists("token.rds")) {
    test_that("resultsTables : Dropbox mode", {
        resTbl <- resultsTables(
            results = res_small,
            counts = dds_small,
            lfcThreshold = lfc,
            summary = FALSE,
            write = TRUE,
            dir = "resultsTables",
            dropboxDir = file.path("bcbioRNASeq_examples", "resultsTables"),
            rdsToken = "token.rds"
        )
        expect_true("dropboxFiles" %in% names(resTbl))
        # Check for Dropbox URLs
        expect_true(all(vapply(
            X = resTbl[["dropboxFiles"]],
            FUN = function(file) {
                grepl("^https://www.dropbox.com/s/", file[["url"]])
            },
            FUN.VALUE = logical(1L)
        )))
        # Now check the Markdown code
        output <- capture.output(.markdownResultsTables(resTbl))
        # The function currently returns links to 4 DEG files
        expect_identical(
            length(which(grepl("https://www.dropbox.com/s/", output))),
            4L
        )
    })
}
