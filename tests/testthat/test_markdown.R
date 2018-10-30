context("Markdown")

data(deseq, envir = environment())
res_small <- as(deseq, "DESeqResults")
res_tables <- DESeqResultsTables(res_small)



# prepareRNASeqTemplate ========================================================
test_that("prepareRNASeqTemplate", {
    files <- c(
        "_footer.Rmd",
        "_header.Rmd",
        "_output.yaml",
        "_setup.R",
        "bibliography.bib"
    )
    expect_silent(prepareRNASeqTemplate())
    expect_true(all(file.exists(files)))
    unlink(files)
})



# topTables ====================================================================
with_parameters_test_that(
    "topTables", {
        x <- capture.output(topTables(object))
        expect_true(grepl("padj", x[[3L]]))
    },
    object = list(
        DESeqAnalysis = deseq,
        DESeqResults = res_small,
        DESeqResultsTables = res_tables
    )
)
