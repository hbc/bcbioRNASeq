data(bcb, envir = environment())



# bcb ====================================================================
context("bcbioRNASeq : bcb")

object <- bcb

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
    },
    object = as.list(assays(object)),
    # nolint start
    sum = list(
        counts = 143758911,
        tpm = 5249973,
        avgTxLength = 2696963,
        normalized = 137992370,
        vst = 32371,
        fpkm = 3800855
    )
    # nolint end
)

rm(object)



# Generator ====================================================================
context("bcbioRNASeq : Generator")

uploadDir <- system.file("extdata/bcbio", package = "bcbioRNASeq")
organism <- "Mus musculus"
ensemblRelease <- 90L

test_that("bcbioRNASeq : salmon (default)", {
    # Expecting warnings about rowRanges.
    object <- suppressWarnings(bcbioRNASeq(uploadDir))
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "tpm", "avgTxLength", "normalized", "vst")
    )

    # Dimensions.
    expect_identical(
        object = dim(object),
        expected = c(100L, 6L)
    )
    expect_identical(
        object = colnames(object),
        expected = paste0(
            rep(c("control_rep", "fa_day7_rep"), each = 3L),
            rep(seq_len(3L), times = 2L)
        )
    )

    # Check that the assays contain expected counts.
    expect_identical(
        object = round(colSums(counts(object))),
        expected = c(
            # nolint start
            control_rep1 = 17659,
            control_rep2 = 60245,
            control_rep3 = 16105,
            fa_day7_rep1 = 29482,
            fa_day7_rep2 = 26272,
            fa_day7_rep3 = 29883
            # nolint end
        )
    )

    # Check that empty ranges were stashed.
    expect_identical(
        object = ncol(rowData(object)),
        expected = 0L
    )
    expect_identical(
        object = levels(seqnames(object)),
        expected = "unknown"
    )

    # Detecting organism automatically if possible.
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = organism
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
            control_rep1 = 21629,
            control_rep2 = 64024,
            control_rep3 = 18773,
            fa_day7_rep1 = 30890,
            fa_day7_rep2 = 28140,
            fa_day7_rep3 = 31869
            # nolint end
        )
    )
})

# Testing both gene and transcript level.
with_parameters_test_that(
    "bcbioRNASeq : GTF/GFF file", {
        skip_on_travis()
        # GFF3 files are also supported, but we're only testing GTF here for
        # speed. This functionality is covered in basejump tests also.
        gtfURL <- paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-90",
            "gtf",
            "mus_musculus",
            "Mus_musculus.GRCm38.90.gtf.gz",
            sep = "/"
        )
        gtfFile <- basename(gtfURL)
        if (!file.exists(gtfFile)) {
            download.file(url = gtfURL, destfile = gtfFile)
        }
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            organism = organism,
            gffFile = gtfFile
        )
        expect_s4_class(object, "bcbioRNASeq")
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
        expected = c(
            "counts",
            "tpm",
            "avgTxLength",
            "normalized"
        )
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
        expected = c(
            "counts",
            "tpm",
            "avgTxLength",
            "normalized",
            "vst",
            "rlog"
        )
    )
})

test_that("bcbioRNASeq : User-defined sample metadata", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        sampleMetadataFile = file.path(uploadDir, "sample_metadata.csv")
    )
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(dim(object), c(100L, 6L))
    expect_identical(
        object = basename(metadata(object)[["sampleMetadataFile"]]),
        expected = "sample_metadata.csv"
    )
})

test_that("bcbioRNASeq : Sample selection", {
    keep <- paste0("control_rep", seq_len(3L))
    censor <- paste0("fa_day7_rep", seq_len(3L))

    # Keep only control samples.
    object <- bcbioRNASeq(uploadDir = uploadDir, samples = keep)
    expect_identical(
        object = colnames(object),
        expected = keep
    )

    # Censor the folic acid samples.
    object <- bcbioRNASeq(uploadDir, censorSamples = censor)
    expect_identical(
        object = colnames(object),
        expected = keep
    )
})



# Methods ======================================================================
context("bcbioRNASeq : Methods")

object <- bcb
nrow <- 50L
ncol <- 2L

test_that("extract : Normal gene and sample selection", {
    subset <- object[seq_len(nrow), seq_len(ncol)]
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = dim(subset),
        expected = c(nrow, ncol)
    )
    expect_identical(
        object = rownames(subset),
        expected = head(rownames(object), n = nrow)
    )
    expect_identical(
        object = colnames(subset),
        expected = head(colnames(object), n = ncol)
    )

    # Require at least 50 genes, 2 samples.
    expect_error(object[seq_len(49L), ])
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
            .[seq_len(nrow), ] %>%
            assayNames(),
        expected = c(
            "counts",
            "tpm",
            "avgTxLength",
            "normalized",
            "vst",
            "fpkm"
        )
    )

    # Allow the user to skip, using `recalculate` argument.
    expect_identical(
        object = object %>%
            .[seq_len(nrow), , recalculate = FALSE] %>%
            assayNames(),
        expected = c(
            "counts",
            "tpm",
            "avgTxLength"
        )
    )
})



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



test_that("updateObject", {
    # Load a legacy object that doesn't contain rowRanges.
    load("bcb_invalid.rda")
    object <- bcb_invalid

    expect_error(
        object = validObject(object),
        regexp = "rowRanges"
    )
    expect_identical(
        object = slot(object, "metadata")[["version"]],
        expected = package_version("0.1.4")
    )

    x <- suppressWarnings(updateObject(object))
    expect_s4_class(x, "bcbioRNASeq")
    expect_true(validObject(x))
    expect_identical(
        object = metadata(x)[["version"]],
        expected = packageVersion
    )
    expect_identical(
        object = metadata(x)[["previousVersion"]],
        expected = package_version("0.1.4")
    )
})
