context("loadRNASeq")

# Load the minimal example bcbio run saved in the package
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
bcb <- loadRNASeq(uploadDir)
annotable <- annotable(bcb)

test_that("Class definition", {
    expect_identical(
        slotNames(bcb),
        c("bcbio",
          "colData",
          "assays",
          "NAMES",
          "elementMetadata",
          "metadata")
    )
    expect_identical(
        slot(bcb, "bcbio") %>% class(),
        structure("SimpleList", package = "S4Vectors")
    )
    expect_identical(
        slot(bcb, "colData") %>% class(),
        structure("DataFrame", package = "S4Vectors")
    )
    expect_identical(
        slot(bcb, "assays") %>% class(),
        structure("ShallowSimpleListAssays", package = "SummarizedExperiment")
    )
    expect_identical(
        slot(bcb, "NAMES") %>% class(),
        "character"
    )
    expect_identical(
        slot(bcb, "elementMetadata") %>% class(),
        structure("DataFrame", package = "S4Vectors")
    )
    # Switch this structure to SimpleList?
    expect_identical(
        slot(bcb, "metadata") %>% class(),
        "list"
    )
})

test_that("Assays", {
    expect_identical(
        lapply(assays(bcb), class),
        list(raw = "matrix",
             normalized = "matrix",
             tpm = "matrix",
             tmm = "matrix",
             rlog = structure("DESeqTransform", package = "DESeq2"),
             vst = structure("DESeqTransform", package = "DESeq2")
        )
    )
})

test_that("Column data", {
    expect_identical(
        lapply(colData(bcb), class),
        list(sampleID = "factor",
             sampleName = "factor",
             description = "factor",
             group = "factor")
    )
})

# Ensembl annotations from AnnotationHub, using ensembldb
test_that("Row data", {
    expect_identical(
        lapply(rowData(bcb), class),
        list(ensgene = "character",
             symbol = "character",
             description = "character",
             biotype = "character",
             broadClass = "character",
             geneSeqStart = "integer",
             geneSeqEnd = "integer",
             seqName = "character",
             seqStrand = "integer",
             seqCoordSystem = "character",
             entrez = "list")
    )
})

test_that("Metadata", {
    expect_identical(
        lapply(metadata(bcb), class),
        list(version = c("package_version", "numeric_version"),
             uploadDir = "character",
             sampleDirs = "character",
             projectDir = "character",
             template = "character",
             runDate = "Date",
             interestingGroups = "character",
             organism = "character",
             genomeBuild = "character",
             ensemblVersion = "NULL",
             annotable = "data.frame",
             tx2gene = "data.frame",
             lanes = "integer",
             yaml = "list",
             metrics = "data.frame",
             sampleMetadataFile = "NULL",
             dataVersions = c("tbl_df", "tbl", "data.frame"),
             programVersions = c("tbl_df", "tbl", "data.frame"),
             bcbioLog = "character",
             bcbioCommandsLog = "character",
             allSamples = "logical",
             design = "formula",
             transformationLimit = "integer",
             date = "Date",
             wd = "character",
             utilsSessionInfo = "sessionInfo",
             devtoolsSessionInfo = "session_info",
             unannotatedGenes = "character")
    )
    # Interesting groups should default to `sampleName`
    expect_identical(
        metadata(bcb)[["interestingGroups"]],
        "sampleName"
    )
    # Ensembl metadata version should default to `NULL`
    expect_identical(
        metadata(bcb)[["ensemblVersion"]],
        NULL
    )
})

test_that("bcbio", {
    expect_identical(
        lapply(slot(bcb, "bcbio"), class),
        list(tximport = "list",
             DESeqDataSet = structure("DESeqDataSet", package = "DESeq2"),
             featureCounts = "matrix")
    )
})

test_that("Example data dimensions", {
    expect_identical(
        dim(bcb),
        c(505L, 4L)
    )
    expect_identical(
        colnames(bcb),
        c("group1_1", "group1_2", "group2_1", "group2_2")
    )
    expect_identical(
        rownames(bcb)[1L:4L],
        c("ENSMUSG00000002459",
          "ENSMUSG00000004768",
          "ENSMUSG00000005886",
          "ENSMUSG00000016918")
    )
})

test_that("transformationLimit", {
    # TODO Remove warning suppression when example dataset is updated
    skip <- suppressWarnings(
        loadRNASeq(
            uploadDir = uploadDir,
            annotable = annotable,
            transformationLimit = 0L)
    )
    expect_identical(
        names(assays(skip)),
        c("raw", "normalized", "tpm", "tmm")
    )
    expect_identical(metadata(skip)[["transformationLimit"]], 0L)
    expect_warning(
            loadRNASeq(
                uploadDir = uploadDir,
                annotable = annotable,
                transformationLimit = 0L),
        paste(
            "Dataset contains many samples.",
            "Skipping DESeq2 variance stabilization."
        )
    )
    expect_message(
        suppressWarnings(
            loadRNASeq(
                uploadDir = uploadDir,
                annotable = annotable,
                transformationLimit = Inf)
        ),
        "Performing rlog transformation"
    )
})

test_that("Metadata from the cloud", {
    sampleMetadataFile <- "http://bcbiornaseq.seq.cloud/sample_metadata.csv"
    bcb <- loadRNASeq(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile)
    expect_s4_class(bcb, "bcbioRNASeq")
    expect_identical(
        metadata(bcb)[["sampleMetadataFile"]],
        sampleMetadataFile
    )
})
