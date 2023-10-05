test_that("Import salmon counts (default)", {
    object <- bcbioRNASeq(uploadDir)
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
    expect_identical(
        object = assayNames(object),
        expected = c(
            "counts",
            "aligned",
            "avgTxLength",
            "tpm",
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

test_that("Fast mode in R", {
    object <- bcbioRNASeq(uploadDir = uploadDir, fast = TRUE)
    expect_s4_class(object, "bcbioRNASeq")
    expect_true(validObject(object))
    expect_identical(
        object = assayNames(object),
        expected = c(
            "counts",
            "avgTxLength",
            "tpm"
        )
    )
})

test_that("bcbio fastrnaseq pipeline", {
    uploadDir <- file.path(cacheDir, "fastrnaseq")
    object <- bcbioRNASeq(uploadDir = uploadDir, fast = TRUE)
    expect_s4_class(object, "bcbioRNASeq")
})

test_that("countsFromAbundance", {
    colSums <- c(
        control_rep1 = 17658.87,
        control_rep2 = 60245.14,
        control_rep3 = 16104.51,
        fa_day7_rep1 = 29481.62,
        fa_day7_rep2 = 26272.19,
        fa_day7_rep3 = 29883.19
    )
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "genes",
        countsFromAbundance = "lengthScaledTPM"
    )
    expect_true(isSubset("avgTxLength", assayNames(x)))
    expect_identical(
        object = round(colSums(assay(x)), digits = 2L),
        expected = colSums
    )
    expect_identical(
        object = round(assay(x)[1L, , drop = TRUE], digits = 2L),
        expected = c(
            control_rep1 = 473.71,
            control_rep2 = 1860.57,
            control_rep3 = 458.73,
            fa_day7_rep1 = 1106.95,
            fa_day7_rep2 = 1048.79,
            fa_day7_rep3 = 1089.97
        )
    )
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "genes",
        countsFromAbundance = "no"
    )
    expect_true(isSubset("avgTxLength", assayNames(x)))
    expect_identical(
        object = round(colSums(assay(x)), digits = 2L),
        expected = colSums
    )
    expect_identical(
        object = round(assay(x)[1L, , drop = TRUE], digits = 2L),
        expected = c(
            ## nolint start
            control_rep1 = 472,
            control_rep2 = 1832,
            control_rep3 = 473,
            fa_day7_rep1 = 1094,
            fa_day7_rep2 = 1051,
            fa_day7_rep3 = 1092
            ## nolint end
        )
    )
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "transcripts",
        countsFromAbundance = "lengthScaledTPM"
    )
    expect_true("avgTxLength" %in% assayNames(x))
    expect_identical(
        object = round(colSums(assay(x)), digits = 2L),
        expected = colSums
    )
    expect_identical(
        object = round(assay(x)[1L, , drop = TRUE], digits = 2L),
        expected = c(
            control_rep1 = 473.75,
            control_rep2 = 1860.89,
            control_rep3 = 459.31,
            fa_day7_rep1 = 1110.09,
            fa_day7_rep2 = 1051.12,
            fa_day7_rep3 = 1093.84
        )
    )
    x <- bcbioRNASeq(
        uploadDir = uploadDir,
        level = "transcripts",
        countsFromAbundance = "no"
    )
    expect_true("avgTxLength" %in% assayNames(x))
    expect_identical(
        object = round(colSums(assay(x)), digits = 2L),
        expected = colSums
    )
    expect_identical(
        object = round(colSums(assay(x)), digits = 2L),
        expected = c(
            control_rep1 = 17658.87,
            control_rep2 = 60245.14,
            control_rep3 = 16104.51,
            fa_day7_rep1 = 29481.62,
            fa_day7_rep2 = 26272.19,
            fa_day7_rep3 = 29883.19
        )
    )
})

## Testing both gene and transcript level.
test_that("AnnotationHub", {
    for (level in eval(formals(bcbioRNASeq)[["level"]])) {
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            caller = "salmon",
            organism = organism,
            ensemblRelease = ensemblRelease
        )
        expect_s4_class(object, "bcbioRNASeq")
    }
})

## Aligned counts can only be loaded at gene level.
test_that("STAR aligned counts", {
    object <- bcbioRNASeq(
        uploadDir = uploadDir,
        caller = "star",
        level = "genes"
    )
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

## Testing both gene and transcript level. GFF3 files are also supported, but
## we're only testing GTF here for speed. This functionality is covered in
## basejump tests also.
test_that("GTF/GFF file", {
    for (level in eval(formals(bcbioRNASeq)[["level"]])) {
        gffURL <- paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-90",
            "gtf",
            "mus_musculus",
            "Mus_musculus.GRCm38.90.gtf.gz",
            sep = "/"
        )
        gffFile <- file.path(cacheDir, basename(gffURL))
        if (!file.exists(gffFile)) {
            initDir(cacheDir)
            download.file(url = gffURL, destfile = gffFile)
        }
        object <- bcbioRNASeq(
            uploadDir = uploadDir,
            level = level,
            organism = organism,
            gffFile = gffFile
        )
        expect_s4_class(object, "bcbioRNASeq")
    }
})
