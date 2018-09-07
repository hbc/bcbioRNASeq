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
    object <- aggregateCols(object)
    expect_identical(
        object = colnames(object),
        expected = c("group1", "group2")
    )
    expect_identical(
        object = sum(counts(object)),
        expected = sum(counts(bcb))
    )
    expect_equal(
        object = rowSums(counts(object)),
        expected = rowSums(counts(bcb))
    )
})



# counts =======================================================================
with_parameters_test_that(
    "counts : Slotted assays", {
        # Check that all are matrices.
        expect_is(
            object = counts(bcb, normalized = normalized),
            class = "matrix"
        )
        # Check that we're matching the expected assay matrix.
        expect_identical(
            object = counts(bcb, normalized = normalized),
            expected = assays(bcb)[[assay]]
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
        expect_is(
            object = counts(bcb, normalized = normalized),
            class = "matrix"
        )
        expect_null(assays(bcb)[[normalized]])
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



# interestingGroups ============================================================
test_that("interestingGroups : NULL validity checks", {
    object <- bcb
    expect_error(interestingGroups(object) <- NULL)
    metadata(object)[["interestingGroups"]] <- NULL
    expect_error(validObject(object))
})



# sampleData ===================================================================
test_that("sampleData : Verbose mode (default)", {
    object <- sampleData(bcb, clean = FALSE)

    # Return `interestingGroups` factor column by default.
    expect_is(object[["interestingGroups"]], "factor")

    # Otherwise it should be identical to `colData`.
    object[["interestingGroups"]] <- NULL
    expected <- colData(bcb)
    expect_identical(object, expected)
})

test_that("sampleData : Clean mode", {
    object <- sampleData(bcb, clean = TRUE)
    # Return only useful factor columns.
    expect_identical(
        object = colnames(object),
        expected = c("sampleName", "group")
    )
    # Require that all clean columns are factor.
    invisible(lapply(object, function(object) {
        expect_is(object, "factor")
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
