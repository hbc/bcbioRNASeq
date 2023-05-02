## nolint start
DESeqAnalysis <- DESeqAnalysis::DESeqAnalysis
cacheURL <- pipette::cacheURL
degPatterns <- DEGreport::degPatterns
degPlot <- DEGreport::degPlot
import <- pipette::import
pasteURL <- AcidBase::pasteURL
plotVolcano <- AcidGenerics::plotVolcano
realpath <- AcidBase::realpath
render <- rmarkdown::render
results <- DESeq2::results
resultsNames <- AcidGenerics::resultsNames
resultsTables <- AcidGenerics::resultsTables
significants <- DEGreport::significants
tempdir2 <- AcidBase::tempdir2
## nolint end

## Manuscript uses `loadRemoteData()` instead.
object <- import(
    con = cacheURL(
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
object <- updateObject(object, verbose = TRUE)
bcb <- object
dds <- as.DESeqDataSet(object)
design(dds) <- ~day
dds <- DESeq(dds)
vst <- varianceStabilizingTransformation(dds)
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
alphaThreshold <- 0.05
lfcThreshold <- 0L
resListUnshrunken <- Map(
    f = DESeq2::results,
    contrast = contrasts,
    MoreArgs = list(
        "object" = dds,
        "alpha" = alphaThreshold,
        "lfcThreshold" = lfcThreshold
    )
)
res <- resListUnshrunken[[1L]]
resListShrunken <- Map(
    f = DESeq2::lfcShrink,
    res = resListUnshrunken,
    contrast = contrasts,
    MoreArgs = list(
        "dds" = dds,
        "type" = "normal",
        "alpha" = alphaThreshold,
        "lfcThreshold" = lfcThreshold
    )
)
deseq <- DESeqAnalysis(
    data = dds,
    transform = vst,
    results = resListUnshrunken,
    lfcShrink = resListShrunken
)

test_that("counts", {
    for (normalized in list(FALSE, TRUE, "tpm", "rlog", "vst")) {
        x <- counts(object, normalized = normalized)
        expect_type(x, "double")
    }
})

test_that("saveData and writeCounts", {
    tempdir <- tempdir2()
    unlink(tempdir)
    raw <- counts(object, normalized = FALSE)
    normalized <- counts(object, normalized = TRUE)
    tpm <- counts(object, normalized = "tpm")
    rlog <- counts(object, normalized = "rlog")
    vst <- counts(object, normalized = "vst")
    saveData(
        raw, normalized, tpm, rlog, vst,
        dir = file.path(tempdir, "saveData")
    )
    expect_identical(
        object = sort(list.files(
            path = file.path(tempdir, "saveData"),
            full.names = FALSE
        )),
        expected = c(
            "normalized.rds",
            "raw.rds",
            "rlog.rds",
            "tpm.rds",
            "vst.rds"
        )
    )
    writeCounts(
        raw, normalized, tpm, rlog, vst,
        dir = file.path(tempdir, "writeCounts")
    )
    expect_identical(
        object = sort(list.files(
            path = file.path(tempdir, "writeCounts"),
            full.names = FALSE
        )),
        expected = c(
            "normalized.csv",
            "raw.csv",
            "rlog.csv",
            "tpm.csv",
            "vst.csv"
        )
    )
    unlink(tempdir)
})

test_that("Quality control", {
    for (fun in list(
        plotTotalReads,
        plotMappingRate,
        plotExonicMappingRate,
        plotIntronicMappingRate,
        plotGenesDetected,
        plotGeneSaturation,
        plotCountsPerGene,
        plotCountsPerFeature,
        plotCountDensity,
        plotMeanSD,
        plotPCA
    )) {
        p <- fun(object)
        expect_s3_class(p, "ggplot")
    }
    p <- plotCorrelationHeatmap(object)
    expect_s3_class(p, "pheatmap")
    for (fun in list(
        plotDispEsts,
        plotPCACovariates
    )) {
        p <- fun(object)
        expect_type(p, "list")
    }
})

test_that("plotCountsPerFeature", {
    p <- plotCountsPerFeature(object, geom = "density")
    expect_identical(nrow(p[["data"]]), 457704L)
    expect_match(p[["labels"]][["subtitle"]], "38142")
})

test_that("Differential expression", {
    x <- capture.output({
        alphaSummary(
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
    ## NOTE The generic formals have been renamed here.
    ## object : results
    ## DESeqTransform : counts
    ## Still providing legacy, deprecated support.
    p <- plotDEGHeatmap(
        object = res,
        DESeqTransform = vst
    )
    expect_s3_class(p, "pheatmap")
    p <- degPlot(
        object,
        res = res,
        n = 3L,
        slot = "vst",
        log2 = FALSE,
        ann = c("geneId", "geneName"),
        xs = "day"
    )
    expect_s3_class(p, "ggplot")
})

test_that("Detecting patterns", {
    ddsLrt <- DESeq(dds, test = "LRT", reduced = ~1L)
    resLrt <- results(ddsLrt)
    ma <- counts(object, "vst")[significants(res, fc = 2L), ]
    resPatterns <- degPatterns(
        ma = ma,
        metadata = colData(object),
        time = "day",
        minc = 60L
    )
    p <- resPatterns[["plot"]]
    expect_s3_class(p, "ggplot")
})

test_that("resultsTables and topTables", {
    ## "lfc" argument renamed to "lfcThreshold".
    resTbl <- resultsTables(res, lfcThreshold = 1L)
    expect_s4_class(resTbl, "DataFrameList")
    expect_s4_class(resTbl[[1L]], "DESeqResults")
    expect_output(topTables(resTbl, n = 5L))
})

## Modified, updated version of bcbio_rnaseq_output_example repo.
## https://github.com/bcbio/bcbio_rnaseq_output_example/

templatesDir <- system.file(
    "rmarkdown",
    "templates",
    package = .pkgName,
    mustWork = TRUE
)
renderDir <- tempdir2()
renderFiles <- saveData(
    bcb, dds, deseq, res,
    dir = renderDir,
    ext = "rds",
    overwrite = TRUE
)
withr::with_dir(
    new = renderDir,
    code = {
        AcidMarkdown::prepareTemplate()
    }
)

test_that("Quality Control", {
    stem <- "01-quality-control"
    input <- file.path(renderDir, paste0(stem, ".Rmd"))
    file.copy(
        from = file.path(
            templatesDir,
            stem,
            "skeleton",
            "skeleton.Rmd"
        ),
        to = input,
        overwrite = TRUE
    )
    x <- render(
        input = input,
        params = list(
            "bcb_file" = renderFiles[["bcb"]]
        ),
        clean = TRUE
    )
    outfile <- file.path(renderDir, paste0(stem, ".html"))
    expect_true(file.exists(outfile))
    expect_identical(realpath(x), realpath(outfile))
})

test_that("Differential Expression", {
    stem <- "02-differential-expression"
    input <- file.path(renderDir, paste0(stem, ".Rmd"))
    file.copy(
        from = file.path(
            templatesDir,
            stem,
            "skeleton",
            "skeleton.Rmd"
        ),
        to = input,
        overwrite = TRUE
    )
    x <- render(
        input = input,
        params = list(
            "bcb_file" = renderFiles[["bcb"]],
            "design" = design(dds),
            "contrasts" = contrasts
        ),
        clean = TRUE
    )
    outfile <- file.path(renderDir, paste0(stem, ".html"))
    expect_true(file.exists(outfile))
    expect_identical(realpath(x), realpath(outfile))
})

## NOTE This test may fail due to ggrepel overlap warning.
## # ggrepel: 5 unlabeled data points (too many overlaps).
## # Consider increasing max.overlaps

test_that("Functional Analysis", {
    stem <- "03-functional-analysis"
    input <- file.path(renderDir, paste0(stem, ".Rmd"))
    file.copy(
        from = file.path(
            templatesDir,
            stem,
            "skeleton",
            "skeleton.Rmd"
        ),
        to = input,
        overwrite = TRUE
    )
    ## Enabling pathview here is too slow for CI.
    x <- render(
        input = input,
        params = list(
            "deseq_file" = renderFiles[["deseq"]],
            "results_name" = resultsNames(deseq)[[1L]],
            "pathview" = FALSE
        )
    )
    outfile <- file.path(renderDir, paste0(stem, ".html"))
    expect_identical(x, outfile)
    expect_true(file.exists(outfile))
})

AcidBase::unlink2(renderDir)
