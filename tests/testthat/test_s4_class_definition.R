context("S4 Class Definition")

bcb <- suppressWarnings(loadRNASeq(
    uploadDir = uploadDir,
    ensemblRelease = 87L,
    organism = "Mus musculus"
))
validObject(bcb)



# bcbioRNASeq S4 Object ========================================================
test_that("Slots", {
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

test_that("Dimensions", {
    expect_identical(
        dim(bcb),
        c(502L, 4L)
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
            "entrezID" = "list",
            "broadClass" = "factor"
        )
    )
})

test_that("Metadata", {
    tibble <- c("tbl_df", "tbl", "data.frame")
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
            "rowRangesMetadata" = tibble,
            "gffFile" = "character",
            "tx2gene" = "data.frame",
            "lanes" = "integer",
            "yaml" = "list",
            "metrics" = "data.frame",
            "dataVersions" = tibble,
            "programVersions" = tibble,
            "bcbioLog" = "character",
            "bcbioCommandsLog" = "character",
            "allSamples" = "logical",
            "loadRNASeq" = "call",
            "date" = "Date",
            "wd" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info"
        )
    )
    # Interesting groups should default to `sampleName`
    expect_identical(
        metadata(bcb)[["interestingGroups"]],
        "sampleName"
    )
})



# subset =======================================================================
test_that("subset : Normal gene and sample selection", {
    x <- bcb_small[seq_len(100L), seq_len(4L)]
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(dim(x), c(100L, 4L))
    expect_identical(
        rownames(x)[[1L]],
        rownames(bcb_small)[[1L]]
    )
    expect_identical(
        colnames(x),
        head(colnames(bcb_small), 4L)
    )
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized", "rlog", "vst")
    )
})

test_that("subset : Minimal selection ranges", {
    # Require at least 100 genes, 2 samples
    x <- bcb_small[seq_len(100L), seq_len(2L)]
    expect_error(bcb_small[seq_len(99L), ])
    expect_error(bcb_small[, seq_len(1L)])
    expect_identical(
        dimnames(x),
        list(
            head(rownames(bcb_small), 100L),
            head(colnames(bcb_small), 2L)
        )
    )
})



# loadRNASeq ===================================================================
test_that("loadRNASeq : organism = NULL", {
    x <- loadRNASeq(
        uploadDir = uploadDir,
        organism = NULL
    )
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        levels(seqnames(x)),
        "unknown"
    )
})

test_that("loadRNASeq : transformationLimit", {
    x <- suppressWarnings(
        loadRNASeq(
            uploadDir = uploadDir,
            organism = "Mus musculus",
            transformationLimit = -Inf
        )
    )
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized")
    )
})

test_that("loadRNASeq : User-defined sample metadata", {
    x <- suppressWarnings(loadRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        sampleMetadataFile = file.path(uploadDir, "sample_metadata.csv")
    ))
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        basename(metadata(x)[["sampleMetadataFile"]]),
        "sample_metadata.csv"
    )
})



# updateObject =================================================================
test_that("updateObject", {
    expect_error(validObject(bcb_invalid))
    expect_identical(
        slot(bcb_invalid, "metadata")[["version"]],
        package_version("0.1.4")
    )
    organism <- slot(bcb_invalid, "metadata")[["organism"]]
    rowRanges <- makeGRangesFromEnsembl(organism)
    x <- suppressWarnings(
        updateObject(bcb_invalid, rowRanges = rowRanges)
    )
    expect_identical(
        metadata(x)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(x)[["previousVersion"]],
        package_version("0.1.4")
    )
})
