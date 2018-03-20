context("loadRNASeq")

bcb <- suppressWarnings(loadRNASeq(
    uploadDir,
    ensemblRelease = 87L,
    organism = "Mus musculus"
))
validObject(bcb)

test_that("Class definition", {
    expect_identical(
        slotNames(bcb),
        c(
            "rowRanges",
            "colData",
            "assays",
            "NAMES",
            "elementMetadata",
            "metadata"
        )
    )
    expect_identical(
        lapply(seq_along(slotNames(bcb)), function(a) {
            class(slot(bcb, slotNames(bcb)[[a]]))
        }),
        list(
            structure(
                "GRanges",
                package = "GenomicRanges"
            ),
            structure(
                "DataFrame",
                package = "S4Vectors"
            ),
            structure(
                "ShallowSimpleListAssays",
                package = "SummarizedExperiment"
            ),
            "NULL",  # character for SummarizedExperiment
            structure(
                "DataFrame",
                package = "S4Vectors"
            ),
            "list"
        )
    )
})

test_that("Assays", {
    # All assays should be a matrix
    expect_true(all(vapply(
        X = assays(bcb),
        FUN = function(assay) {
            is.matrix(assay)
        },
        FUN.VALUE = logical(1L)
    )))
})

test_that("Column data", {
    # All columns should be factor5
    expect_true(all(vapply(
        X = colData(bcb),
        FUN = function(assay) {
            is.factor(assay)
        },
        FUN.VALUE = logical(1L)
    )))
})

# TODO Add rowRanges check

test_that("Row data", {
    # Ensembl annotations from AnnotationHub, using ensembldb
    expect_identical(
        lapply(rowData(bcb), class),
        list(
            "geneID" = "character",
            "geneName" = "character",
            "geneBiotype" = "factor",
            "description" = "character",
            "seqCoordSystem" = "factor",
            "entrezID" = "AsIs",
            "broadClass" = "factor"
        )
    )
})

test_that("Metadata", {
    tbl_df <- c("tbl_df", "tbl", "data.frame")
    expect_identical(
        lapply(metadata(bcb), class),
        list(
            "version" = c("package_version", "numeric_version"),
            "level" = "character",
            "caller" = "character",
            "countsFromAbundance" = "character",
            "uploadDir" = "character",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "projectDir" = "character",
            "template" = "character",
            "runDate" = "Date",
            "interestingGroups" = "character",
            "organism" = "character",
            "genomeBuild" = "character",
            "ensemblRelease" = "integer",
            "rowRangesMetadata" = tbl_df,
            "gffFile" = "character",
            "tx2gene" = "data.frame",
            "lanes" = "integer",
            "yaml" = "list",
            "metrics" = "data.frame",
            "dataVersions" = tbl_df,
            "programVersions" = tbl_df,
            "bcbioLog" = "character",
            "bcbioCommandsLog" = "character",
            "allSamples" = "logical",
            "loadRNASeq" = "call",
            "date" = "Date",
            "wd" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info",
            "isSpike" = "character",
            "unannotatedRows" = "character"
        )
    )
    # Interesting groups should default to `sampleName`
    expect_identical(
        metadata(bcb)[["interestingGroups"]],
        "sampleName"
    )
})

test_that("Example data dimensions", {
    expect_identical(
        dim(bcb),
        c(503L, 4L)
    )
    expect_identical(
        colnames(bcb),
        c("group1_1", "group1_2", "group2_1", "group2_2")
    )
    expect_identical(
        rownames(bcb)[1L:4L],
        c(
            "ENSMUSG00000002459",
            "ENSMUSG00000004768",
            "ENSMUSG00000005886",
            "ENSMUSG00000016918"
        )
    )
})

test_that("transformationLimit", {
    skip <- suppressWarnings(
        loadRNASeq(
            uploadDir = uploadDir,
            organism = "Mus musculus",
            transformationLimit = -Inf
        )
    )
    expect_identical(
        names(assays(skip)),
        c("raw", "tpm", "length")
    )
})

test_that("User-defined sample metadata", {
    bcb <- suppressWarnings(loadRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        sampleMetadataFile = sampleMetadataFile
    ))
    expect_s4_class(bcb, "bcbioRNASeq")
    expect_identical(
        metadata(bcb)[["sampleMetadataFile"]],
        sampleMetadataFile
    )
})
