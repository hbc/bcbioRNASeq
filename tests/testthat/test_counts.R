context("counts")

load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))

test_that("raw", {
    counts <- counts(bcb, normalized = FALSE)
    expect_identical(
        assay(bcb),
        counts
    )
})

test_that("normalized", {
    counts <- counts(bcb, normalized = TRUE)
    expect_identical(
        assays(bcb)[["normalized"]],
        counts
    )
})

test_that("tpm", {
    counts <- counts(bcb, normalized = "tpm")
    tpm <- tpm(bcb)
    expect_identical(counts, tpm)
})

test_that("tmm", {
    counts <- counts(bcb, normalized = "tmm")
    tmm <- tmm(bcb)
    expect_identical(counts, tmm)
})

test_that("rlog", {
    counts <- counts(bcb, normalized = "rlog")
    dt <- assays(bcb)[["rlog"]]
    rlog <- assay(dt)
    expect_identical(rlog, counts)
})

test_that("vst", {
    counts <- counts(bcb, normalized = "vst")
    dt <- assays(bcb)[["vst"]]
    vst <- assay(dt)
    expect_identical(vst, counts)
})

test_that("transformationLimit", {
    skip <- bcb
    assays(skip)[["rlog"]] <- NULL
    expect_warning(
        counts(skip, normalized = "rlog"),
        paste(
            "rlog counts not defined.",
            "Calculating and using log2 tmm counts on the fly instead."
        )
    )
    counts <- suppressWarnings(counts(skip, normalized = "rlog"))
    expect_is(counts, "matrix")
})
