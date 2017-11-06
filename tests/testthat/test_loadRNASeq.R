context("loadRNASeq")

# Load the minimal example bcbio run saved in the package
uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
# This should produce the following warnings:
#   1: bcbio-nextgen.log missing
#   2: bcbio-nextgen-commands.log missing
#   3: Unannotated genes detected in counts matrix
bcb <- suppressWarnings(loadRNASeq(uploadDir))

test_that("class definition", {
    expect_equal(
        slotNames(bcb),
        c("bcbio",
          "colData",
          "assays",
          "NAMES",
          "elementMetadata",
          "metadata")
    )
    expect_equal(
        slot(bcb, "bcbio") %>% class(),
        structure("SimpleList", package = "S4Vectors")
    )
    expect_equal(
        slot(bcb, "colData") %>% class(),
        structure("DataFrame", package = "S4Vectors")
    )
    expect_equal(
        slot(bcb, "assays") %>% class(),
        structure("ShallowSimpleListAssays",
                  package = "SummarizedExperiment")
    )
    expect_equal(
        slot(bcb, "NAMES") %>% class(),
        "character"
    )
    expect_equal(
        slot(bcb, "elementMetadata") %>% class(),
        structure("DataFrame", package = "S4Vectors")
    )
    # Switch this structure to SimpleList?
    expect_equal(
        slot(bcb, "metadata") %>% class(),
        "list"
    )
})

test_that("assays", {
    expect_equal(
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

test_that("colData", {
    expect_equal(
        lapply(colData(bcb), class),
        list(sampleID = "factor",
             sampleName = "factor",
             description = "factor",
             group = "factor")
    )
})

# Ensembl annotations from AnnotationHub, using ensembldb
test_that("rowData", {
    expect_equal(
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
             seqCoordSystem = "character")
    )
})

test_that("metadata", {
    expect_equal(
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
             lanes = "numeric",
             yaml = "list",
             metrics = "data.frame",
             sampleMetadataFile = "NULL",
             dataVersions = c("tbl_df", "tbl", "data.frame"),
             programs = c("tbl_df", "tbl", "data.frame"),
             bcbioLog = "NULL",
             bcbioCommandsLog = "NULL",
             allSamples = "logical",
             date = "Date",
             wd = "character",
             utilsSessionInfo = "sessionInfo",
             devtoolsSessionInfo = "session_info",
             unannotatedGenes = "character")
    )
    # Interesting groups should default to `sampleName`
    expect_equal(
        metadata(bcb)[["interestingGroups"]],
        "sampleName"
    )
    # Ensembl metadata version should default to `NULL`
    expect_equal(
        metadata(bcb)[["ensemblVersion"]],
        NULL
    )
})

test_that("bcbio", {
    expect_equal(
        lapply(slot(bcb, "bcbio"), class),
        list(tximport = "list",
             DESeqDataSet = structure("DESeqDataSet", package = "DESeq2"),
             featureCounts = "matrix")
    )
})

test_that("example data dimensions", {
    expect_equal(
        dim(bcb),
        c(505, 4)
    )
    expect_equal(
        colnames(bcb),
        c("group1_1", "group1_2", "group2_1", "group2_2")
    )
    expect_equal(
        rownames(bcb)[1:4],
        c("ENSMUSG00000002459",
          "ENSMUSG00000004768",
          "ENSMUSG00000005886",
          "ENSMUSG00000016918")
    )
})
