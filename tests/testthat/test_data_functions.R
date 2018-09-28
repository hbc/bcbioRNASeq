context("Data Functions")

bcb <- bcbioRNASeq(
    uploadDir = system.file("extdata/bcbio", package = "bcbioRNASeq"),
    rlog = TRUE,
    vst = TRUE
)
dds <- as(bcb, "DESeqDataSet")



# aggregateCols ================================================================
test_that("aggregateCols", {
    object <- bcb
    # Assign groupings into `aggregate` column of `colData()`.
    aggregate <- as.factor(sub("^([a-z0-9]+)_.*", "\\1", colnames(object)))
    names(aggregate) <- colnames(object)
    object[["aggregate"]] <- aggregate
    x <- aggregateCols(object)
    expect_identical(
        object = colnames(x),
        expected = c("group1", "group2")
    )
    expect_identical(
        object = sum(counts(x)),
        expected = sum(counts(object))
    )
    expect_equal(
        object = rowSums(counts(x)),
        expected = rowSums(counts(object))
    )
})



# counts =======================================================================
with_parameters_test_that(
    "counts : Slotted assays", {
        object <- bcb
        # Check that all are matrices.
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        # Check that we're matching the expected assay matrix.
        expect_identical(
            object = counts(object, normalized = normalized),
            expected = assays(object)[[assay]]
        )
    },
    normalized = list(
        FALSE,
        TRUE,
        "tpm",
        "vst",
        "rlog"
    ),
    assay = list(
        "counts",
        "normalized",
        "tpm",
        "vst",
        "rlog"
    )
)

with_parameters_test_that(
    "counts : On the fly assays", {
        object <- bcb
        expect_is(
            object = counts(object, normalized = normalized),
            class = "matrix"
        )
        expect_null(assays(object)[[normalized]])
    },
    normalized = c("tmm", "rle")
)

with_parameters_test_that(
    "counts : Skipped DESeq2 transforms", {
        object <- bcb
        assays(object)[[normalized]] <- NULL
        expect_warning(
            object = counts(object, normalized = normalized),
            regexp = paste(
                "Calculating log2 TMM counts instead"
            )
        )
        expect_is(
            object = suppressWarnings(
                counts(object, normalized = normalized)
            ),
            class = "matrix"
        )
    },
    normalized = c("rlog", "vst")
)



# sampleData ===================================================================
test_that("sampleData : Verbose mode (default)", {
    object <- bcb
    x <- sampleData(object, clean = FALSE)

    # Return `interestingGroups` factor column by default.
    expect_is(x[["interestingGroups"]], "factor")

    # Otherwise it should be identical to `colData`.
    x[["interestingGroups"]] <- NULL
    expect_identical(
        object = x,
        expected = colData(object)
    )
})

test_that("sampleData : Clean mode", {
    object <- bcb
    x <- sampleData(object, clean = TRUE)
    # Return only useful factor columns.
    expect_identical(
        object = colnames(x),
        expected = c("sampleName", "group")
    )
    # Require that all clean columns are factor.
    invisible(lapply(x, function(x) {
        expect_is(x, "factor")
    }))
})

test_that("sampleData<-", {
    object <- bcb
    sampleData(object)[["testthat"]] <- factor("XXX")
    expect_identical(
        object = levels(sampleData(object)[["testthat"]]),
        expected = "XXX"
    )
})



# selectSamples ================================================================
with_parameters_test_that(
    "selectSamples", {
        # All objects have the same `group` column defined.
        object <- selectSamples(object, group = "ctrl")
        expect_s4_class(object, class = "SummarizedExperiment")
        expect_identical(
            object = colnames(object),
            expected = c("group1_1", "group1_2")
        )
    },
    object = list(
        bcbioRNASeq = bcb,
        DESeqDataSet = dds
    )
)



# tmm ==========================================================================
with_parameters_test_that(
    "tmm", {
        expect_is(
            object = tmm(object),
            class = "matrix"
        )
    },
    object = list(
        bcbioRNASeq = bcb,
        DESeqDataSet = dds,
        matrix = assay(bcb)
    )
)
