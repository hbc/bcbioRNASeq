context("bcbioRNASeq")

test_that("Import salmon counts (default)", {
    object <- bcbioRNASeq(uploadDir)
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
    expect_identical(
        object = assayNames(object),
        expected = c(
            "counts",
            "tpm",
            "avgTxLength",
            "aligned",
            "normalized",
            "vst"
        )
    )

    ## Dimensions.
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

    ## Check that the assays contain expected counts.
    expect_identical(
        object = round(colSums(counts(object))),
        expected = c(
            ## nolint start
            control_rep1 = 17659,
            control_rep2 = 60245,
            control_rep3 = 16105,
            fa_day7_rep1 = 29482,
            fa_day7_rep2 = 26272,
            fa_day7_rep3 = 29883
            ## nolint end
        )
    )

    ## Check that empty ranges were stashed.
    expect_identical(
        object = ncol(rowData(object)),
        expected = 0L
    )
    expect_identical(
        object = levels(seqnames(object)),
        expected = "unknown"
    )

    ## Detecting organism automatically if possible.
    expect_identical(
        object = metadata(object)[["organism"]],
        expected = organism
    )
})

test_that("Fast mode", {
    object <- bcbioRNASeq(uploadDir = uploadDir, fast = TRUE)
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
    expect_identical(
        object = assayNames(object),
        expected = c(
            "counts",
            "tpm",
            "avgTxLength"
        )
    )
})

## Testing both gene and transcript level.
with_parameters_test_that(
    "AnnotationHub", {
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

## Aligned counts can only be loaded at gene level.
test_that("STAR aligned counts", {
    object <- bcbioRNASeq(uploadDir, caller = "star")
    expect_s4_class(object, "bcbioRNASeq")
    expect_identical(
        object = assayNames(object),
        expected = c("counts", "normalized", "vst")
    )
    ## Aligned counts are integer.
    expect_identical(
        object = colSums(counts(object)),
        expected = c(
            ## nolint start
            control_rep1 = 21629,
            control_rep2 = 64024,
            control_rep3 = 18773,
            fa_day7_rep1 = 30890,
            fa_day7_rep2 = 28140,
            fa_day7_rep3 = 31869
            ## nolint end
        )
    )
})

test_that("User-defined sample metadata", {
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

test_that("Sample selection", {
    keep <- paste0("control_rep", seq_len(3L))
    censor <- paste0("fa_day7_rep", seq_len(3L))

    ## Keep only control samples.
    object <- bcbioRNASeq(uploadDir = uploadDir, samples = keep)
    expect_identical(
        object = colnames(object),
        expected = keep
    )

    ## Censor the folic acid samples.
    object <- bcbioRNASeq(uploadDir, censorSamples = censor)
    expect_identical(
        object = colnames(object),
        expected = keep
    )
})

## Testing both gene and transcript level.
with_parameters_test_that(
    "GTF/GFF file", {
        skip_on_appveyor()
        skip_on_docker()
        ## GFF3 files are also supported, but we're only testing GTF here for
        ## speed. This functionality is covered in basejump tests also.
        gtfURL <- paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-90",
            "gtf",
            "mus_musculus",
            "Mus_musculus.GRCm38.90.gtf.gz",
            sep = "/"
        )
        gtfFile <- file.path("cache", basename(gtfURL))
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
