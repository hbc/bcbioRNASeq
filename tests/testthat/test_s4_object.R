context("S4 Object")

bcb <- bcbioRNASeq(
    uploadDir = uploadDir,
    organism = "Mus musculus",
    ensemblRelease = 87L
)



# bcbioRNASeq Object Structure =================================================
test_that("Slot names", {
    expect_identical(
        object = slotNames(bcb),
        expected = c(
            "rowRanges",
            "colData",
            "assays",
            "NAMES",
            "elementMetadata",
            "metadata"
        )
    )
})

with_parameters_test_that(
    "Slot definitions", {
        expect_identical(
            object = class(slot(bcb, slotName)),
            expected = expected
        )
    },
    slotName = slotNames(bcb),
    expected = list(
        rowRanges = structure(
            "GRanges",
            package = "GenomicRanges"
        ),
        colData = structure(
            "DataFrame",
            package = "S4Vectors"
        ),
        assays = structure(
            "ShallowSimpleListAssays",
            package = "SummarizedExperiment"
        ),
        NAMES = "NULL",
        elementMetadata = structure(
            "DataFrame",
            package = "S4Vectors"
        ),
        metadata = "list"
    )
)

test_that("Dimensions", {
    expect_identical(
        object = dim(bcb),
        expected = c(502L, 4L)
    )
    expect_identical(
        object = colnames(bcb),
        expected = c("group1_1", "group1_2", "group2_1", "group2_2")
    )
    expect_identical(
        object = head(rownames(bcb), n = 4L),
        expected = c(
            "ENSMUSG00000002459",
            "ENSMUSG00000004768",
            "ENSMUSG00000005886",
            "ENSMUSG00000016918"
        )
    )
})

with_parameters_test_that(
    "Assays", {
        expect_is(object, "matrix")
        # Check to make sure tximport loaded correctly.
        # Transcript-to-gene counts aren't integer, so we're rounding here to
        # check the values more easily.
        expect_identical(
            object = round(sum(object)),
            expected = sum
        )
        expect_identical(
            object = round(colSums(object)),
            expected = colSums
        )
    },
    object = as.list(assays(bcb)),
    # nolint start
    sum = list(
        counts = 1861378,
        tpm = 145541,
        length = 2587914,
        normalized = 1837956,
        vst = 13069
    ),
    colSums = list(
        counts = c(
            group1_1 = 289453,
            group1_2 = 515082,
            group2_1 = 494089,
            group2_2 = 562753
        ),
        tpm = c(
            group1_1 = 26800,
            group1_2 = 36631,
            group2_1 = 38619,
            group2_2 = 43491
        ),
        length = c(
            group1_1 = 646685,
            group1_2 = 638827,
            group2_1 = 656783,
            group2_2 = 645620
        ),
        normalized = c(
            group1_1 = 319030,
            group1_2 = 499308,
            group2_1 = 472147,
            group2_2 = 547471
        ),
        vst = c(
            group1_1 = 3266,
            group1_2 = 3266,
            group2_1 = 3272,
            group2_2 = 3265
        )
    )
    # nolint end
)

# Ensembl annotations from AnnotationHub, using ensembldb.
with_parameters_test_that(
    "Row data structure", {
        expect_identical(object, expected)
    },
    object = lapply(rowData(bcb), class),
    expected = list(
        broadClass = "factor",
        description = "factor",
        entrezID = "list",
        geneBiotype = "factor",
        geneID = "character",
        geneName = "factor",
        seqCoordSystem = "factor"
    )
)

test_that("Metadata", {
    tibble <- c("tbl_df", "tbl", "data.frame")
    expect_identical(
        object = lapply(metadata(bcb), class),
        expected = list(
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
            rowRangesMetadata = tibble,
            gffFile = "character",
            tx2gene = structure("tx2gene", package = "basejump"),
            lanes = "integer",
            yaml = "list",
            dataVersions = tibble,
            programVersions = tibble,
            bcbioLog = "character",
            bcbioCommandsLog = "character",
            allSamples = "logical",
            call = "call",
            date = "Date",
            wd = "character",
            utilsSessionInfo = "sessionInfo",
            devtoolsSessionInfo = "session_info"
        )
    )
})

test_that("Metadata values", {
    # Interesting groups should default to `sampleName`.
    expect_identical(
        object = metadata(bcb)[["interestingGroups"]],
        expected = "sampleName"
    )
})



# bcbioRNASeq Constructor ======================================================
test_that("bcbioRNASeq : Aligned counts", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        caller = "star"
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "normalized", "vst")
    )
    # Aligned counts are integer.
    expect_identical(
        object = sum(counts(object)),
        expected = 979020L
    )
    expect_identical(
        object = colSums(counts(object)),
        expected = c(
            # nolint start
            group1_1 = 222940,
            group1_2 = 255196,
            group2_1 = 250838,
            group2_2 = 250046
            # nolint end
        )
    )
})

test_that("bcbioRNASeq : Transcripts", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "transcripts"
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "tpm", "length")
    )
    # Transcript-level counts are not integer.
    # FIXME Ensure we're returning with the rownames sorted.
})

test_that("bcbioRNASeq : organism = NULL", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        organism = NULL
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        levels(seqnames(object)),
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
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        gffFile = gtfFile
    )
})

test_that("bcbioRNASeq : DESeq2 variance stabilization", {
    object <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = FALSE,
            rlog = FALSE
        )
    )
    expect_identical(
        names(assays(object)),
        c("counts", "tpm", "length", "normalized")
    )
    object <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = TRUE,
            rlog = TRUE
        )
    )
    expect_identical(
        names(assays(object)),
        c("counts", "tpm", "length", "normalized", "vst", "rlog")
    )
})

test_that("bcbioRNASeq : User-defined sample metadata", {
    object <- suppressWarnings(bcbioRNASeq(
        uploadDir = uploadDir,
        organism = "Mus musculus",
        sampleMetadataFile = file.path(uploadDir, "sample_metadata.csv")
    ))
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        basename(metadata(object)[["sampleMetadataFile"]]),
        "sample_metadata.csv"
    )
})

test_that("bcbioRNASeq: Sample selection", {
    # samples
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        samples = head(colnames(bcb), 2L)
    )
    expect_identical(
        colnames(object),
        head(colnames(bcb), 2L)
    )

    # censor samples
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        censorSamples = head(colnames(bcb), 1L)
    )
    expect_identical(
        colnames(object),
        colnames(bcb)[-1L]
    )
})



# extract ======================================================================
test_that("extract : Normal gene and sample selection", {
    object <- bcb[seq_len(100L), seq_len(4L)]
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(dim(object), c(100L, 4L))
    expect_identical(
        rownames(object)[[1L]],
        rownames(bcb)[[1L]]
    )
    expect_identical(
        colnames(object),
        head(colnames(bcb), 4L)
    )
    expect_identical(
        assayNames(object),
        c("counts", "tpm", "length", "normalized", "vst")
    )
})

test_that("extract : Minimal selection ranges", {
    # Require at least 100 genes, 2 samples
    object <- bcb[seq_len(100L), seq_len(2L)]
    expect_error(bcb[seq_len(99L), ])
    expect_error(bcb[, seq_len(1L)])
    expect_identical(
        dimnames(object),
        list(
            head(rownames(bcb), 100L),
            head(colnames(bcb), 2L)
        )
    )
})

test_that("extract : DESeq2 transforms", {
    # Transform by default
    object <- bcb[1L:100L, 1L:2L]
    expect_identical(
        assayNames(object),
        c("counts", "tpm", "length", "normalized", "vst")
    )

    # Allow the user to skip, using `transform` argument
    object <- bcb[1L:100L, 1L:2L, transform = FALSE]
    expect_identical(
        names(assays(object)),
        c("counts", "tpm", "length", "normalized")
    )
})

test_that("extract : unmodified", {
    object <- bcb[, ]
    expect_identical(object, bcb)
})



# show =========================================================================
test_that("show", {
    object <- capture.output(show(bcb))
    expect_true(grepl("bcbioRNASeq", object[[1L]]))

    # Fake metadata for code coverage
    object <- bcb
    metadata(object)[["sampleMetadataFile"]] <- "XXX"
    metadata(object)[["gffFile"]] <- "XXX"
    object <- capture.output(show(object))
    expect_true(grepl("bcbioRNASeq", object[[1L]]))
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
    object <- suppressWarnings(updateObject(bcb_invalid))
    expect_s4_class(object, "bcbioRNASeq")
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
    object <- suppressWarnings(
        updateObject(bcb_invalid, rowRanges = rowRanges)
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        metadata(object)[["version"]],
        packageVersion
    )
    expect_identical(
        metadata(object)[["previousVersion"]],
        package_version("0.1.4")
    )
})
