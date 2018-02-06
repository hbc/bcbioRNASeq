context("resultsTables")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioRNASeq"))
load(system.file(
    file.path("extdata", "res.rda"),
    package = "bcbioRNASeq"))

test_that("resultsTables", {
    resTbl <- suppressMessages(
        resultsTables(
        res,
        lfc = 0.25,
        summary = FALSE,
        write = FALSE)
    )
    expect_identical(class(resTbl), "list")
    tibble <- c("tbl_df", "tbl", "data.frame")
    expect_identical(
        lapply(resTbl, class),
        list("contrast" = "character",
             "alpha" = "numeric",
             "lfc" = "numeric",
             "all" = tibble,
             "deg" = tibble,
             "degLFC" = tibble,
             "degLFCUp" = tibble,
             "degLFCDown" = tibble)
    )
})

test_that("Summary and write support", {
    annotable <- annotable(bcb)
    # This is also capturing the contents of the list return. Not sure how to
    # fix here for knitr asis_output
    output <- capture.output(
        resultsTables(
        res,
        lfc = 0.25,
        annotable = annotable,
        summary = TRUE,
        headerLevel = 2L,
        write = TRUE,
        dir = "resultsTables"),
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
            "- 505 genes in counts matrix",
            "- Base mean > 0: 304 genes (non-zero)",
            "- Base mean > 1: 213 genes",
            "- Alpha cutoff: 0.1",
            "- LFC cutoff: 0.25 (applied in tables only)",
            "- DEG pass alpha: 4 genes",
            "- DEG LFC up: 2 genes",
            "- DEG LFC down: 2 genes",
            "",
            "",
            "",
            "## Results tables",
            "",
            "",
            paste0(
                "- [`group_ko_vs_ctrl_all.csv.gz`]",
                "(resultsTables/group_ko_vs_ctrl_all.csv.gz): ",
                "All genes, sorted by Ensembl identifier."
            ),
            paste0(
                "- [`group_ko_vs_ctrl_deg.csv.gz`]",
                "(resultsTables/group_ko_vs_ctrl_deg.csv.gz): ",
                "Genes that pass the alpha (FDR) cutoff."
            ),
            paste0(
                "- [`group_ko_vs_ctrl_deg_lfc_up.csv.gz`]",
                "(resultsTables/group_ko_vs_ctrl_deg_lfc_up.csv.gz): ",
                "Upregulated DEG; positive log2 fold change."
            ),
            paste0(
                "- [`group_ko_vs_ctrl_deg_lfc_down.csv.gz`]",
                "(resultsTables/group_ko_vs_ctrl_deg_lfc_down.csv.gz): ",
                "Downregulated DEG; negative log2 fold change."
            ),
            ""
        )
    )
    # Check for the output files
    files <- list.files("resultsTables")
    expect_identical(
        files,
        c(
            "group_ko_vs_ctrl_all.csv.gz",
            "group_ko_vs_ctrl_deg_lfc_down.csv.gz",
            "group_ko_vs_ctrl_deg_lfc_up.csv.gz",
            "group_ko_vs_ctrl_deg.csv.gz"
        )
    )
    unlink("resultsTables", recursive = TRUE)
})
