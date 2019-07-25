context("updateObject")

skip_if_not(hasInternet())

## Load a legacy object that doesn't contain rowRanges.
load(file.path("cache", "bcb_invalid.rda"))

test_that("bcb_invalid", {
    expect_error(
        object = validObject(bcb_invalid),
        regexp = "rowRanges"
    )
    expect_identical(
        object = slot(bcb_invalid, "metadata")[["version"]],
        expected = package_version("0.1.4")
    )
})

test_that("Expected success", {
    x <- updateObject(bcb_invalid)
    expect_s4_class(x, "bcbioRNASeq")
    expect_true(validObject(x))
    expect_identical(
        object = metadata(x)[["version"]],
        expected = .version
    )
    expect_identical(
        object = metadata(x)[["previousVersion"]],
        expected = package_version("0.1.4")
    )
})

test_that("metadata slot updates", {
    ## genomeBuild
    metadata(bcb_invalid)[["genomeBuild"]] <- FALSE

    ## gtf
    metadata(bcb_invalid)[["gtf"]] <- TRUE

    ## gtfFile
    gffFile <- "XXX.gtf.gz"
    metadata(bcb_invalid)[["gtfFile"]] <- gffFile

    ## missingGenes
    missingGenes <- "XXX"

    ## yamlFile
    yamlFile <- "XXX"

    x <- updateObject(bcb_invalid)
    expect_true(validObject(x))

    ## genomeBuild
    expect_identical(
        object = metadata(x)[["genomeBuild"]],
        expected = character()
    )

    ## gtf
    expect_null(metadata(x)[["gtf"]])

    ## gtfFile
    expect_null(metadata(x)[["gtfFile"]])
    expect_identical(metadata(x)[["gffFile"]], gffFile)

    ## missingGenes
    expect_null(metadata(x)[["missingGenes"]])

    ## yamlFile
    expect_null(metadata(x)[["yamlFile"]])
})
