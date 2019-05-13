context("S4 Object")

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")



# bcbioRNASeq Object Structure =================================================
bcb <- suppressWarnings(bcbioRNASeq(
    uploadDir = uploadDir,
    organism = "Mus musculus",
    ensemblRelease = 87L,
    testthat = "XXX"
))
expect_true(validObject(bcb))

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

test_that("Row data", {
    # Ensembl annotations from AnnotationHub, using ensembldb
    expect_identical(
        lapply(rowData(bcb), class),
        list(
            broadClass = "factor",
            description = "factor",
            entrezID = "list",
            geneBiotype = "factor",
            geneID = "character",
            geneName = "factor",
            seqCoordSystem = "factor"
        )
    )
})

test_that("Metadata", {
    expect_identical(
        lapply(metadata(bcb), class),
        list(
            version = c("package_version", "numeric_version"),
            level = "character",
            caller = "character",
            countsFromAbundance = "character",
            uploadDir = "character",
            sampleDirs = "character",
            sampleMetadataFile = "character",
            projectDir = "character",
            template = "character",
            runDate = "Date",
            interestingGroups = "character",
            organism = "character",
            genomeBuild = "character",
            ensemblRelease = "integer",
            rowRangesMetadata = c("tbl_df", "tbl", "data.frame"),
            gffFile = "character",
            tx2gene = "data.frame",
            lanes = "integer",
            yaml = "list",
            dataVersions = c("tbl_df", "tbl", "data.frame"),
            # `spec_tbl_df` is a new virtual class in tibble v2.0 update.
            programVersions = c("spec_tbl_df", "tbl_df", "tbl", "data.frame"),
            bcbioLog = "character",
            bcbioCommandsLog = "character",
            allSamples = "logical",
            call = "call",
            # Check for user-defined metadata
            testthat = "character",
            date = "Date",
            wd = "character",
            utilsSessionInfo = "sessionInfo",
            devtoolsSessionInfo = "session_info"
        )
    )
    # Interesting groups should default to `sampleName`
    expect_identical(
        metadata(bcb)[["interestingGroups"]],
        "sampleName"
    )
})



# bcbioRNASeq ==================================================================
test_that("bcbioRNASeq : Aligned counts", {
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        caller = "star"
    )
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        assayNames(x),
        c("counts", "normalized", "vst")
    )
})

test_that("bcbioRNASeq : Transcripts", {
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "transcripts"
    )
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        assayNames(x),
        c("counts", "tpm", "length")
    )
})

test_that("bcbioRNASeq : organism = NULL", {
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        organism = NULL
    )
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        levels(seqnames(x)),
        "unknown"
    )
})

# GFF3 files are also supported, but we're not testing for speed
test_that("bcbioRNASeq : GTF file", {
    gtfURL <- paste(
        "ftp://ftp.ensembl.org",
        "pub",
        "release-87",
        "gtf",
        "mus_musculus",
        "Mus_musculus.GRCm38.87.gtf.gz",
        sep = "/"
    )
    gtfFile <- basename(gtfURL)
    if (!file.exists(gtfFile)) {
        download.file(url = gtfURL, destfile = gtfFile)
    }
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        gffFile = gtfFile
    )
})

test_that("bcbioRNASeq : DESeq2 variance stabilization", {
    x <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = FALSE,
            rlog = FALSE
        )
    )
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized")
    )
    x <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = TRUE,
            rlog = TRUE
        )
    )
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized", "vst", "rlog")
    )
})

test_that("bcbioRNASeq : User-defined sample metadata", {
    x <- suppressWarnings(bcbioRNASeq(
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

test_that("bcbioRNASeq: Sample selection", {
    # samples
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        samples = head(colnames(bcb), 2L)
    )
    expect_identical(
        colnames(x),
        head(colnames(bcb), 2L)
    )

    # censor samples
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        censorSamples = head(colnames(bcb), 1L)
    )
    expect_identical(
        colnames(x),
        colnames(bcb)[-1L]
    )
})



# extract ======================================================================
test_that("extract : Normal gene and sample selection", {
    x <- bcb[seq_len(100L), seq_len(4L)]
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(dim(x), c(100L, 4L))
    expect_identical(
        rownames(x)[[1L]],
        rownames(bcb)[[1L]]
    )
    expect_identical(
        colnames(x),
        head(colnames(bcb), 4L)
    )
    expect_identical(
        assayNames(x),
        c("counts", "tpm", "length", "normalized", "vst")
    )
})

test_that("extract : Minimal selection ranges", {
    # Require at least 100 genes, 2 samples
    x <- bcb[seq_len(100L), seq_len(2L)]
    expect_error(bcb[seq_len(99L), ])
    expect_error(bcb[, seq_len(1L)])
    expect_identical(
        dimnames(x),
        list(
            head(rownames(bcb), 100L),
            head(colnames(bcb), 2L)
        )
    )
})

test_that("extract : DESeq2 transforms", {
    # Transform by default
    x <- bcb[1L:100L, 1L:2L]
    expect_identical(
        assayNames(x),
        c("counts", "tpm", "length", "normalized", "vst")
    )

    # Allow the user to skip, using `transform` argument
    x <- bcb[1L:100L, 1L:2L, transform = FALSE]
    expect_identical(
        names(assays(x)),
        c("counts", "tpm", "length", "normalized")
    )
})

test_that("extract : unmodified", {
    x <- bcb[, ]
    expect_identical(x, bcb)
})



# show =========================================================================
test_that("show", {
    x <- capture.output(show(bcb))
    expect_true(grepl("bcbioRNASeq", x[[1L]]))

    # Fake metadata for code coverage
    object <- bcb
    metadata(object)[["sampleMetadataFile"]] <- "XXX"
    metadata(object)[["gffFile"]] <- "XXX"
    x <- capture.output(show(object))
    expect_true(grepl("bcbioRNASeq", x[[1L]]))
})



# updateObject =================================================================
test_that("updateObject", {
    load("bcb_invalid.rda")
    expect_error(validObject(bcb_invalid))
    expect_identical(
        slot(bcb_invalid, "metadata")[["version"]],
        package_version("0.1.4")
    )

    # NULL rowRanges (default)
    x <- suppressWarnings(updateObject(bcb_invalid))
    expect_s4_class(x, "bcbioRNASeq")
    expect_warning(
        updateObject(bcb_invalid),
        "`rowRanges` are now recommended for gene annotations"
    )

    # Rich rowRanges metadata
    organism <- slot(bcb_invalid, "metadata")[["organism"]]
    expect_is(organism, "character")
    rowRanges <- makeGRangesFromEnsembl(organism, release = 87L)
    # Suppressing expected warning about genes:
    # ENSMUSG00000104475, ENSMUSG00000109048
    x <- suppressWarnings(
        updateObject(bcb_invalid, rowRanges = rowRanges)
    )
    expect_s4_class(x, "bcbioRNASeq")
    expect_identical(
        metadata(x)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(x)[["previousVersion"]],
        package_version("0.1.4")
    )
})
