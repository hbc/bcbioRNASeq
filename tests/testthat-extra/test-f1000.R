suppressPackageStartupMessages({
    library(rmarkdown)
    library(DESeq2)
    library(DESeqAnalysis)
})



context("F1000 workflow paper")

## Note that DESeqAnalysis approach is new and not mentioned in the current
## F1000v2 manuscript.

## bcbioRNASeq object.
bcb <- import(
    file = cacheURL(
        url = pasteURL(
            "github.com",
            "hbc",
            "bcbioRNASeq",
            "raw",
            "f1000v2/data/bcb.rda",
            protocol = "https"
        ),
        pkg = .pkgName
    )
)
bcb <- updateObject(bcb)
expect_s4_class(bcb, "bcbioRNASeq")
## Note that "object" is used instead of "bcb" in some tests below.
object <- bcb

alphaThreshold <- 0.05
lfcThreshold <- 0L

## Run DESeq2.
dds <- as(bcb, "DESeqDataSet")
design(dds) <- ~day
dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)

## > resultsNames(dds)
## [1] "Intercept"  "day_1_vs_0" "day_3_vs_0" "day_7_vs_0"

contrasts <- list(
    "day_1_vs_0" = c(
        "factor" = "day",
        "numerator" = 1L,
        "denominator" = 0L
    ),
    "day_3_vs_0" = c(
        "factor" = "day",
        "numerator" = 3L,
        "denominator" = 0L
    ),
    "day_7_vs_0" = c(
        "factor" = "day",
        "numerator" = 7L,
        "denominator" = 0L
    )
)

resListUnshrunken <- mapply(
    FUN = DESeq2::results,
    contrast = contrasts,
    MoreArgs = list(
        "object" = dds,
        "alpha" = alphaThreshold,
        "lfcThreshold" = lfcThreshold
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
)
names(resListUnshrunken) <- names(contrasts)

res <- resListUnshrunken[[1L]]

resListShrunken <- mapply(
    FUN = DESeq2::lfcShrink,
    res = resListUnshrunken,
    contrast = contrasts,
    MoreArgs = list(
        "dds" = dds,
        "type" = "normal",
        "alpha" = alphaThreshold,
        "lfcThreshold" = lfcThreshold
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE
)

## See `help(topic = DESeqAnalysis, DESeqAnalysis)` for details.
deseq <- DESeqAnalysis(
    data = dds,                   # DESeqDataSet
    transform = vst,              # DESeqTransform
    results = resListUnshrunken,  # DESeqResults list (unshrunken)
    lfcShrink = resListShrunken   # DESeqResults list (shrunken)
)


## FIXME Need to save this in renderFiles....see below.


test_that("counts", {
    for (normalized in list(FALSE, TRUE, "tpm", "rlog", "vst"
    )) {
        x <- counts(object, normalized = normalized)
        expect_is(x, "matrix")
    }
})

test_that("Quality control", {
    for (fun in list(
        plotCountDensity,
        plotCountsPerFeature,
        plotCountsPerGene,
        plotExonicMappingRate,
        plotGeneSaturation,
        plotGenesDetected,
        plotIntronicMappingRate,
        plotMappingRate,
        plotMeanSD,
        plotPCA,
        plotTotalReads
    )) {
        p <- fun(object)
        expect_s3_class(p, "ggplot")
    }

    p <- plotCorrelationHeatmap(object)
    expect_s3_class(p, "pheatmap")

    p <- plotCountsPerFeature(object, geom = "density")
    expect_identical(nrow(p[["data"]]), 457704L)
    expect_match(p[["labels"]][["subtitle"]], "38142")

    ## NOTE Figure margins can error here inside of testthat call.
    ## This may be removed in a future update, supporting only DESeqDataSet.
    ## > p <- plotDispEsts(object)
    ## > expect_is(p, "list")

    ## NOTE Disabled until bug is fixed in DEGreport.
    ## > p <- plotPCACovariates(object, fdr = 0.1)
    ## > expect_is(p, "list")
})

test_that("Differential expression", {
    x <- capture.output({
        DESeqAnalysis::alphaSummary(
            object = dds,
            contrast = contrasts[["day_7_vs_0"]],
            alpha = c(0.1, 0.05)
        )
    })
    expect_identical(
        object = x[[2L]],
        expected = "LFC > 0 (up)    1521  1102"
    )

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

    ## "lfc" argument renamed to "lfcThreshold".
    resTbl <- resultsTables(res, lfcThreshold = 1L)
    expect_is(resTbl, "list")
    expect_is(resTbl[[1L]], "tbl_df")

    expect_output(topTables(resTbl, n = 5L))
})



context("Render R Markdown templates")

## Modified, updated version of bcbio_rnaseq_output_example repo.
## https://github.com/bcbio/bcbio_rnaseq_output_example/

templatesDir <- system.file(
    "rmarkdown",
    "templates",
    package = .pkgName,
    mustWork = TRUE
)

renderDir <- paste("render", Sys.Date(), sep = "-")
unlink(renderDir, recursive = TRUE)
renderFiles <- saveData(
    bcb, dds, deseq, res,
    dir = renderDir,
    ext = "rds",
    overwrite = TRUE
)
wd <- getwd()
setwd(renderDir)

test_that("Quality Control", {
    file.copy(
        from = file.path(
            templatesDir,
            "01-quality-control",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "quality-control.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "quality-control.Rmd",
        params = list(
            "bcb_file" = renderFiles[["bcb"]]
        ),
        clean = TRUE
    )
    expect_identical(
        object = basename(x),
        expected = "quality-control.html"
    )
})

test_that("Differential Expression", {
    file.copy(
        from = file.path(
            templatesDir,
            "02-differential-expression",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "differential-expression.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "differential-expression.Rmd",
        params = list(
            "bcb_file" = renderFiles[["bcb"]],
            "design" = design(dds),
            "contrasts" = contrasts
        ),
        clean = TRUE
    )
    expect_identical(
        object = basename(x),
        expected = "differential-expression.html"
    )
})

test_that("Functional Analysis", {
    file.copy(
        from = file.path(
            templatesDir,
            "03-functional-analysis",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "functional-analysis.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "functional-analysis.Rmd",
        params = list(
            "deseq_file" = renderFiles[["deseq"]],
            "results_name" = resultsNames(deseq)[[1L]],
            "organism" = organism(as(deseq, "DESeqDataSet")),
            "pathview" = FALSE
        )
    )
    expect_identical(
        object = basename(x),
        expected = "functional-analysis.html"
    )
})

setwd(wd)
