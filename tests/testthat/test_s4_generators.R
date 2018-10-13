# FIXME Rework this and reorganize other unit tests.
# FIXME Check slots using bcb_small for speed



# bcbioRNASeq ==================================================================
context("S4 Generators : bcbioRNASeq")

# FIXME Parameterize and test transcripts also.
# FIXME Parameterize unit test for STAR/featureCounts.

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 87L

# GFF3 files are also supported, but we're only testing GTF here for speed.
# This functionality is provided by basejump and covered by unit tests.
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

test_that("bcbioRNASeq : organism = NULL", {
    # Expecting warnings about rowRanges.
    object <- suppressWarnings(bcbioRNASeq(uploadDir))
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "tpm", "avgTxLength", "normalized", "vst")
    )
    # Detecting organism automatically if possible.
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = organism
    )
    expect_identical(
        object = ncol(rowData(object)),
        expected = 0L
    )
    expect_identical(
        object = levels(seqnames(object)),
        expected = "unknown"
    )
})

# Testing both gene and transcript level.
with_parameters_test_that(
    "bcbioRNASeq : AnnotationHub", {
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            caller = "salmon",
            organism = organism,
            ensemblRelease = ensemblRelease
        )
        expect_s4_class(object, "bcbioRNASeq")
        # tximport counts are length scaled by default and not integer.
        expect_identical(
            object = round(colSums(counts(object))),
            expected = c(
                # nolint start
                group1_1 = 289453,
                group1_2 = 515082,
                group2_1 = 494089,
                group2_2 = 562753
                # nolint end
            )
        )
    },
    level = eval(formals(bcbioRNASeq)[["level"]])
)

# Aligned counts can only be loaded at gene level.
test_that("bcbioRNASeq : Aligned counts", {
    object <- bcbioRNASeq(uploadDir, caller = "star")
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "normalized", "vst")
    )
    # Aligned counts are integer.
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

# Testing both gene and transcript level.
with_parameters_test_that(
    "bcbioRNASeq : GTF/GFF file", {
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            organism = organism,
            gffFile = gtfFile
        )
        expect_s4_class(object, "bcbioRNASeq")
        mcols <- c(
            "broadClass",
            "geneBiotype",
            "geneID",
            "geneName",
            "geneSource",
            "geneVersion",
            "havanaGene",
            "havanaGeneVersion",
            "source",
            "type"
        )
        expect_true(all(mcols %in% colnames(rowData(object))))
    },
    level = eval(formals(bcbioRNASeq)[["level"]])
)

test_that("bcbioRNASeq : DESeq2 variance stabilization", {
    object <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = FALSE,
            rlog = FALSE
        )
    )
    expect_identical(
        object = names(assays(object)),
        expected = c("counts", "tpm", "length", "normalized")
    )

    object <- suppressWarnings(
        bcbioRNASeq(
            uploadDir = uploadDir,
            vst = TRUE,
            rlog = TRUE
        )
    )
    expect_identical(
        object = names(assays(object)),
        expected = c("counts", "tpm", "length", "normalized", "vst", "rlog")
    )
})

test_that("bcbioRNASeq : User-defined sample metadata", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        sampleMetadataFile = file.path(uploadDir, "sample_metadata.csv")
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = basename(metadata(object)[["sampleMetadataFile"]]),
        expected = "sample_metadata.csv"
    )
})

test_that("bcbioRNASeq: Sample selection", {
    keep <- head(colnames(object), n = 2L)
    censor <- setdiff(colnames(object), keep)

    # keep samples
    object <- bcbioRNASeq(uploadDir, samples = keep)
    expect_identical(
        object = colnames(object),
        expected = keep
    )

    # censor samples
    object <- bcbioRNASeq(uploadDir, censorSamples = censor)
    expect_identical(
        object = colnames(object),
        expected = keep
    )
})



# bcbioRNASeq Object Structure =================================================
test_that("Slot names", {
    expect_identical(
        object = slotNames(object),
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
            object = class(slot(object, slotName)),
            expected = expected
        )
    },
    slotName = slotNames(object),
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
        object = dim(object),
        expected = c(502L, 4L)
    )
    expect_identical(
        object = colnames(object),
        expected = c("group1_1", "group1_2", "group2_1", "group2_2")
    )
    expect_identical(
        object = head(rownames(object), n = 4L),
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
    object = as.list(assays(object)),
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
    object = lapply(rowData(object), class),
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
        object = lapply(metadata(object), class),
        expected = list(
            version = c("package_version", "numeric_version"),
            level = "character",
            caller = "character",
            countsFromAbundance = "character",
            uploadDir = "character",
            sampleDirs = "character",
            sampleMetadataFile = "character",
            projectDir = "character",
            runDate = "Date",
            interestingGroups = "character",
            organism = "character",
            genomeBuild = "character",
            ensemblRelease = "integer",
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
            sessionInfo = "session_info"
        )
    )
})

test_that("Metadata values", {
    # Interesting groups should default to `sampleName`.
    expect_identical(
        object = metadata(object)[["interestingGroups"]],
        expected = "sampleName"
    )
})



# extract ======================================================================
test_that("extract : Normal gene and sample selection", {
    subset <- object[seq_len(100L), seq_len(4L)]
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = dim(subset),
        expected = c(100L, 4L)
    )
    expect_identical(
        object = rownames(subset),
        expected = head(rownames(object), n = 100L)
    )
    expect_identical(
        object = colnames(subset),
        expected = head(colnames(object), n = 4L)
    )

    # Require at least 100 genes, 2 samples.
    expect_s4_class(
        object = object[seq_len(100L), seq_len(2L)],
        class = "bcbioRNASeq"
    )
    expect_error(object[seq_len(99L), ])
    expect_error(object[, seq_len(1L)])

    # Check for unmodified return when using empty brackets.
    expect_identical(
        object = object[, ],
        expected = object
    )
})

test_that("extract : Calculate DESeq2 transforms", {
    # Transform enabled by default (if calculated).
    expect_identical(
        object = object %>%
            .[seq_len(100L), ] %>%
            assayNames(),
        expected = c("counts", "tpm", "length", "normalized", "vst")
    )

    # Allow the user to skip, using `recalculate` argument.
    expect_identical(
        object = object %>%
            .[seq_len(100L), , recalculate = FALSE] %>%
            assayNames(),
        expected = c("counts", "tpm", "length")
    )
})



# show =========================================================================
test_that("show", {
    # Stash fake metadata for code coverage.
    metadata(object)[["sampleMetadataFile"]] <- "XXX"
    metadata(object)[["gffFile"]] <- "XXX"
    output <- capture.output(object)
    # Ensure that show method contains "bcbioRNASeq" in the first line.
    expect_true(
        grepl(
            pattern = "^bcbioRNASeq",
            x = output[[1L]]
        )
    )
})



# updateObject =================================================================
test_that("updateObject", {
    # Load a legacy object that doesn't contain rowRanges.
    load("bcb_invalid.rda")
    expect_error(
        object = validObject(bcb_invalid),
        regexp = "rowRanges"
    )
    expect_identical(
        object = slot(bcb_invalid, "metadata")[["version"]],
        expected = package_version("0.1.4")
    )

    # Update using rowRanges.
    organism <- slot(bcb_invalid, "metadata")[["organism"]]
    expect_is(organism, "character")
    rowRanges <- makeGRangesFromEnsembl(organism, release = 87L)
    # Suppressing expected warning about genes:
    # ENSMUSG00000104475, ENSMUSG00000109048
    object <- suppressWarnings(
        updateObject(bcb_invalid, rowRanges = rowRanges)
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
    expect_identical(
        object = metadata(object)[["version"]],
        expected = packageVersion
    )
    expect_identical(
        object = metadata(object)[["previousVersion"]],
        expected = package_version("0.1.4")
    )

    # Update without rowRanges argument (`NULL`).
    # This isn't recommended, but is supported.
    expect_warning(
        object = updateObject(bcb_invalid),
        regexp = "`rowRanges` are now recommended for gene annotations"
    )
    object <- suppressWarnings(updateObject(bcb_invalid))
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
})



# Transcript-Level Counts ======================================================
context("Transcript-Level Counts")

# Loading without rowRanges here, for speed.
object <- bcbioRNASeq(uploadDir, level = "transcripts")

test_that("counts", {
    # Valid: FALSE, tpm.
    expect_is(counts(object, normalized = FALSE), "matrix")
    expect_is(counts(object, normalized = "tpm"), "matrix")
    # All other options are invalid.
    expect_error(counts(object, normalized = TRUE))
    expect_error(counts(object, normalized = "vst"))
})

test_that("Show method", {
    output <- capture.output(show(object))
    expect_true(any(grepl("transcripts", output)))
})

test_that("Invalid functions", {
    expect_error(
        object = plotDispEsts(object),
        regexp = "Gene-level counts are required"
    )
})



# export
# markdown

# FIXME This needs an update.
# test_that("resultsTables : DESeqAnalysis", {
#     output <- capture.output(
#         resultsTables(
#             object = res_small,
#             counts = dds_small,
#             lfcThreshold = lfc,
#             summary = TRUE,
#             headerLevel = 2L,
#             write = TRUE,
#             dir = "resultsTables"
#         ),
#         type = "output"
#     ) %>%
#         # Just evaluate the R Markdown output
#         head(28L)
#     expect_identical(
#         output,
#         c(
#             "",
#             "",
#             "## treatment folic acid vs control",
#             "",
#             "",
#             "",
#             "### Summary statistics",
#             "",
#             "",
#             "- 500 genes in counts matrix",
#             "- Base mean > 0: 500 genes (non-zero)",
#             "- Base mean > 1: 500 genes",
#             "- Alpha: 0.1",
#             "- LFC threshold: 0.25",
#             "- DEG pass alpha: 254 genes",
#             "- DEG LFC up: 115 genes",
#             "- DEG LFC down: 139 genes",
#             "",
#             "",
#             "",
#             "### Results tables",
#             "",
#             "",
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_all.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_all.csv.gz"
#                 ),
#                 "): ",
#                 "All genes, sorted by Ensembl identifier."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg.csv.gz"
#                 ),
#                 "): ",
#                 "Genes that pass the alpha (FDR) cutoff."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg_lfc_up.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg_lfc_up.csv.gz"
#                 ),
#                 "): ",
#                 "Upregulated DEG; positive log2 fold change."
#             ),
#             paste0(
#                 "- [`treatment_folic_acid_vs_control_deg_lfc_down.csv.gz`]",
#                 "(",
#                 file.path(
#                     dir,
#                     "treatment_folic_acid_vs_control_deg_lfc_down.csv.gz"
#                 ),
#                 "): ",
#                 "Downregulated DEG; negative log2 fold change."
#             ),
#             ""
#         )
#     )
# })
#
# if (file.exists("token.rds")) {
#     test_that("resultsTables : DESeqAnalysis : Dropbox mode", {
#         resTbl <- resultsTables(
#             results = res_small,
#             counts = dds_small,
#             lfcThreshold = lfc,
#             summary = FALSE,
#             write = TRUE,
#             dir = "resultsTables",
#             dropboxDir = file.path("bcbioRNASeq_examples", "resultsTables"),
#             rdsToken = "token.rds"
#         )
#         expect_true("dropboxFiles" %in% names(resTbl))
#         # Check for Dropbox URLs
#         expect_true(all(vapply(
#             X = resTbl[["dropboxFiles"]],
#             FUN = function(file) {
#                 grepl("^https://www.dropbox.com/s/", file[["url"]])
#             },
#             FUN.VALUE = logical(1L)
#         )))
#         # Now check the Markdown code
#         output <- capture.output(.markdownResultsTables(resTbl))
#         # The function currently returns links to 4 DEG files
#         expect_identical(
#             length(which(grepl("https://www.dropbox.com/s/", output))),
#             4L
#         )
#     })
# }

