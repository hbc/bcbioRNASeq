context("Differential Expression Functions")

genes <- head(rownames(res_small))
gene2symbol <- gene2symbol(bcb_small)



# alphaSummary =================================================================
test_that("alphaSummary : bcbioRNASeq", {
    expect_warning(
        alphaSummary(bcb_small),
        "Empty design formula detected"
    )
    expect_is(suppressWarnings(alphaSummary(bcb_small)), "knitr_kable")
})

test_that("alphaSummary : DESeqDataSet", {
    x <- alphaSummary(dds_small)
    expect_is(x, "knitr_kable")
    expect_true(grepl("1e-06", x[[1L]]))
})



# plotDEGPCA ===================================================================
test_that("DESeqResults, DESeqTransform", {
    p <- plotDEGPCA(res_small, counts = rld_small)
    expect_is(p, "ggplot")
})



# plotDEGHeatmap ===============================================================
test_that("plotDEGHeatmap", {
    p <- plotDEGHeatmap(res_small, counts = rld_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
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

test_that("plotMA : Gene labels", {
    p <- plotMA(res_small, genes = genes)
    expect_is(p, "ggplot")
})



# plotVolcano ==================================================================
test_that("plotVolcano : DESeqResults", {
    p <- plotVolcano(res_small, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label the top genes
    p <- plotVolcano(res_small, ntop = 5L, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")

    # Label specific genes
    p <- plotVolcano(res_small, genes = genes, gene2symbol = gene2symbol)
    expect_is(p, "ggplot")
})



# resultsTables ================================================================
test_that("resultsTables : Default return with local files only", {
    resTbl <- resultsTables(
        res_small,
        lfc = lfc,
        summary = FALSE,
        write = FALSE
    )
    expect_identical(class(resTbl), "list")
    tibble <- c("tbl_df", "tbl", "data.frame")
    expect_identical(
        lapply(resTbl, class),
        list(
            "contrast" = "character",
             "alpha" = "numeric",
             "lfc" = "numeric",
             "all" = tibble,
             "deg" = tibble,
             "degLFC" = tibble,
             "degLFCUp" = tibble,
             "degLFCDown" = tibble
        )
    )
})

test_that("resultsTables : Summary and write support", {
    # This is also capturing the contents of the list return. Not sure how to
    # fix here for knitr asis_output
    dir <- file.path(getwd(), "resultsTables")
    output <- capture.output(
        resultsTables(
            object = res_small,
            lfc = lfc,
            rowData = rowData(bcb_small),
            summary = TRUE,
            headerLevel = 2L,
            write = TRUE,
            dir = "resultsTables"
        ),
        type = "output"
    ) %>%
        # Just evaluate the R Markdown output
        head(24L)
    expect_identical(
        output,
        c(
            "",
            "",
            "## Summary statistics",
            "",
            "",
            "- 500 genes in counts matrix",
            "- Base mean > 0: 500 genes (non-zero)",
            "- Base mean > 1: 500 genes",
            "- Alpha cutoff: 0.1",
            "- LFC cutoff: 0.25 (applied in tables only)",
            "- DEG pass alpha: 254 genes",
            "- DEG LFC up: 115 genes",
            "- DEG LFC down: 139 genes",
            "",
            "",
            "",
            "## Results tables",
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

if (file.exists("token.rds")) {
    test_that("resultsTables : Dropbox mode", {
        resTbl <- resultsTables(
            res_small,
            lfc = lfc,
            rowData = rowData(bcb_small),
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
