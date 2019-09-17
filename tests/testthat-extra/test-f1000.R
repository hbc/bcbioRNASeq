## nolint start



## Functions.
DESeq <- DESeq2::DESeq
design <- DESeq2::design
`design<-` <- `DESeq2::design<-`
render <- rmarkdown::render
results <- DESeq2::results
varianceStabilizingTransformation <- DESeq2::varianceStabilizingTransformation

## bcbioRNASeq object.
object <- import("https://github.com/hbc/bcbioRNASeq/raw/f1000v2/data/bcb.rda")
object <- updateObject(object)
expect_s4_class(object, "bcbioRNASeq")

## DESeq2 example objects.
bcb <- object
dds <- as(object, "DESeqDataSet")
group <- "day"
contrast <- c(group = group, numerator = "7", denominator = "0")
design_formula <- as.formula(paste0("~", group))
design(dds) <- design_formula
dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
res <- results(object = dds, contrast = contrast, alpha = 0.05)



context("F1000 workflow paper")

test_that("counts", {
    raw <- counts(object, normalized = FALSE)
    expect_is(raw, "matrix")

    normalized <- counts(object, normalized = TRUE)
    expect_is(normalized, "matrix")

    tpm <- counts(object, normalized = "tpm")
    expect_is(tpm, "matrix")

    rlog <- counts(object, normalized = "rlog")
    expect_is(rlog, "matrix")

    vst <- counts(object, normalized = "vst")
    expect_is(vst, "matrix")
})

test_that("Quality control", {
    p <- plotTotalReads(object)
    expect_s3_class(p, "ggplot")

    p <- plotMappingRate(object)
    expect_s3_class(p, "ggplot")

    p <- plotExonicMappingRate(object)
    expect_s3_class(p, "ggplot")

    p <- plotIntronicMappingRate(object)
    expect_s3_class(p, "ggplot")

    expect_warning(
        object = plotGenesDetected(object),
        regexp = "deprecated"
    )

    p <- plotGeneSaturation(object)
    expect_s3_class(p, "ggplot")

    expect_warning(
        object = plotCountsPerGene(object),
        regexp = "deprecated"
    )

    p <- plotCountsPerFeature(object, geom = "density")
    expect_identical(nrow(p$data), 457704L)
    expect_match(p$labels$subtitle, "38142")

    p <- plotCountDensity(object)
    expect_s3_class(p, "ggplot")

    ## Skipping these, they're slow.

    ## > p <- plotMeanSD(object)
    ## > expect_s3_class(p, "ggplot")

    ## > p <- plotDispEsts(object)
    ## > expect_is(p, "list")

    p <- plotCorrelationHeatmap(object)
    expect_s3_class(p, "pheatmap")

    p <- plotPCA(object)
    expect_s3_class(p, "ggplot")

    p <- plotPCACovariates(object, fdr = 0.1)
    expect_is(p, "list")
})

test_that("Differential expression", {
    x <- capture.output(
        alphaSummary(
            object = dds,
            contrast = contrast,
            alpha = c(0.1, 0.05)
        )
    )
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

    ## lfc argument renamed to lfcThreshold.
    res_tbl <- resultsTables(res, lfcThreshold = 1)
    expect_is(res_tbl, "list")
    expect_is(res_tbl[[1L]], "tbl_df")

    expect_output(topTables(res_tbl, n = 5))
})



context("Render R Markdown templates")

## Modified, updated version of bcbio_rnaseq_output_example repo.
## https://github.com/bcbio/bcbio_rnaseq_output_example/

templates_dir <- system.file("rmarkdown", "templates", package = "bcbioRNASeq")
render_dir <- "render"
unlink(render_dir, recursive = TRUE)

render_files <- saveData(
    bcb, dds, res,
    dir = render_dir,
    ext = "rds",
    overwrite = TRUE
)

setwd(render_dir)

test_that("Quality Control", {
    file.copy(
        from = file.path(
            templates_dir,
            "quality-control",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "quality-control.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "quality-control.Rmd",
        params = list(bcb_file = render_files[["bcb"]]),
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
            templates_dir,
            "differential-expression",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "differential-expression.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "differential-expression.Rmd",
        params = list(
            bcb_file = render_files[["bcb"]],
            design = design_formula,
            contrasts = list(contrast)
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
            templates_dir,
            "functional-analysis",
            "skeleton",
            "skeleton.Rmd"
        ),
        to = "functional-analysis.Rmd",
        overwrite = TRUE
    )
    x <- render(
        input = "functional-analysis.Rmd",
        params = list(
            dds_file = dds_file,
            res_file = res_file,
            organism = "Mus musculus",
            data_dir = data_dir
        )
    )
    expect_identical(basename(x), "fa.html")
})

setwd("..")



## nolint end
