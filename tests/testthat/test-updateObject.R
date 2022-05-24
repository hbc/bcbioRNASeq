skip_if_not(hasInternet())

## Load a legacy object that doesn't contain rowRanges.
invalid <- import(file.path(bcbioRNASeqTestsURL, "bcbioRNASeq_0.1.4.rds"))

test_that("v0.1.4 up", {
    expect_error(
        object = validObject(invalid),
        regexp = "rowRanges"
    )
    expect_identical(
        object = slot(invalid, "metadata")[["version"]],
        expected = package_version("0.1.4")
    )
})

test_that("Expected success", {
    x <- updateObject(invalid)
    expect_s4_class(x, "bcbioRNASeq")
    expect_true(validObject(x))
    expect_identical(
        object = metadata(x)[["packageVersion"]],
        expected = .pkgVersion
    )
    expect_identical(
        object = metadata(x)[["previousVersion"]],
        expected = package_version("0.1.4")
    )
})

test_that("metadata slot updates", {
    ## genomeBuild
    metadata(invalid)[["genomeBuild"]] <- FALSE

    ## gtf
    metadata(invalid)[["gtf"]] <- TRUE

    ## gtfFile
    gffFile <- "XXX.gtf.gz"
    metadata(invalid)[["gtfFile"]] <- gffFile

    ## missingGenes
    missingGenes <- "XXX"

    ## yamlFile
    yamlFile <- "XXX"

    x <- updateObject(invalid)
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
