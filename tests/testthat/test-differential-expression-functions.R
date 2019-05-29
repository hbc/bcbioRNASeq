context("Differential Expression Functions")

genes <- head(rownames(res_small))
gene2symbol <- Gene2Symbol(bcb_small)



# alphaSummary =================================================================
test_that("alphaSummary : DESeqDataSet", {
    x <- alphaSummary(dds_small)
    expect_is(x, "knitr_kable")
    expect_true(grepl("1e-06", x[[1L]]))

    x <- alphaSummary(
        object = dds_small,
        contrast = c("treatment", "folic_acid", "control")
    )
    expect_is(x, "knitr_kable")

    x <- alphaSummary(
        object = dds_small,
        name = "treatment_folic_acid_vs_control"
    )
    expect_is(x, "knitr_kable")
})



# plotDEGHeatmap ===============================================================
test_that("plotDEGHeatmap : bcbioRNASeq", {
    p <- plotDEGHeatmap(
        results = res_small,
        counts = bcb_small,
        normalized = "vst"
    )
    expect_identical(names(p), pheatmapNames)
})

test_that("plotDEGHeatmap : DESeqTransform", {
    p <- plotDEGHeatmap(
        results = res_small,
        counts = rld_small
    )
    expect_identical(names(p), pheatmapNames)
})

test_that("plotDEGHeatmap : DESeqDataSet", {
    p <- plotDEGHeatmap(
        results = res_small,
        counts = dds_small,
        title = NULL
    )
    expect_identical(names(p), pheatmapNames)
})

test_that("plotDEGHeatmap : No DEGs", {
    p <- suppressWarnings(plotDEGHeatmap(
        results = res_small,
        counts = bcb_small,
        lfcThreshold = Inf
    ))
    expect_null(p)
})



# plotDEGPCA ===================================================================
test_that("plotDEGPCA : DESeqResults, bcbioRNASeq", {
    p <- plotDEGPCA(results = res_small, counts = bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotDEGPCA : DESeqResults, DESeqTransform", {
    p <- plotDEGPCA(results = res_small, counts = rld_small)
    expect_is(p, "ggplot")
})

test_that("plotDEGPCA : DESeqResults, DESeqDataSet", {
    p <- plotDEGPCA(results = res_small, counts = dds_small)
    expect_is(p, "ggplot")
})

test_that("plotDEGPCA : No DEGs", {
    p <- suppressWarnings(plotDEGPCA(
        results = res_small,
        counts = bcb_small,
        lfcThreshold = Inf
    ))
    expect_null(p)
})



# plotMA =======================================================================
test_that("plotMA : DESeqResults", {
    p <- plotMA(res_small)
    expect_is(p, "ggplot")

    # Check geom classes
    geomtype <- vapply(
        X = p[["layers"]],
        FUN = function(x) {
            class(x[["geom"]])[[1L]]
        },
        FUN.VALUE = character(1L)
    )
    expect_identical(
        geomtype,
        c("GeomHline", "GeomPoint", "GeomLogticks")
    )

    # Check plot labels
    expect_identical(
        p[["labels"]][["y"]],
        "log2 fold change"
    )
    expect_identical(
        p[["labels"]][["x"]],
        "mean expression across all samples"
    )
})

test_that("plotMA : Specific genes", {
    p <- plotMA(
        object = res_small,
        genes = genes,
        gene2symbol = Gene2Symbol(bcb_small)
    )
    expect_is(p, "ggplot")
})

test_that("plotMA: ntop mode", {
    p <- plotMA(
        object = res_small,
        ntop = 10L,
        gene2symbol = Gene2Symbol(bcb_small)
    )
    expect_is(p, "ggplot")
})

test_that("plotMA : Directional support", {
    # Upregulated
    p <- plotMA(
        object = res_small,
        direction = "up",
        sigPointColor = "red"
    )
    expect_is(p, "ggplot")

    # Downregulated
    p <- plotMA(
        object = res_small,
        direction = "down",
        sigPointColor = "green"
    )
    expect_is(p, "ggplot")
})

test_that("plotMA : data.frame return", {
    x <- plotMA(res_small, return = "data.frame")
    expect_is(x, "data.frame")
    expect_true("isDE" %in% colnames(x))
})



# plotVolcano ==================================================================
test_that("plotVolcano : DESeqResults", {
    p <- plotVolcano(res_small, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Enable histograms
    p <- plotVolcano(res_small, histograms = TRUE)
    expect_is(p, "ggplot")

    # Label the top genes
    p <- plotVolcano(res_small, ntop = 5L, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label specific genes
    p <- plotVolcano(res_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Directional support
    p <- plotVolcano(res_small, direction = "up", sigPointColor = "red")
    expect_is(p, "ggplot")
    p <- plotVolcano(res_small, direction = "down", sigPointColor = "green")
    expect_is(p, "ggplot")

    # Return data.frame
    x <- plotVolcano(res_small, return = "data.frame")
    expect_is(x, "data.frame")
})



# resultsTables ================================================================
test_that("resultsTables : Default return with local files only", {
    x <- resultsTables(
        results = res_small,
        counts = dds_small,
        lfcThreshold = lfc,
        summary = FALSE,
        write = FALSE
    )
    expect_identical(class(x), "list")
    expect_identical(
        lapply(x, class),
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
        rownames(x[["all"]])
    )
    expect_identical(
        counts(dds_small, normalized = TRUE),
        as.matrix(x[["all"]][, colnames(assay(dds_small))])
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
    x <- resultsTables(results = res_small)
    expect_is(x, "list")
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
