#' Import Transcript-Level Abundances and Counts
#'
#' Import RNA-seq counts using
#' [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport).
#'
#' Normalized counts are loaded as length-scaled transcripts per million.
#' Consult this [vignette](https://goo.gl/h6fm15) for more information.
#'
#' @note Ignoring transcript versions should work by default. There may be
#' an issue with genomes containing non-Ensembl transcript IDs, such as
#' C. elegans, although we need to double check.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @inheritParams tximport::tximport
#' @param sampleDirs `character`. Sample directories to import.
#' @param type `string`. Expression caller to use.
#'
#' @seealso [tximport::tximport()].
#'
#' @return `list`.
.tximport <- function(
    sampleDirs,
    type = c("salmon", "kallisto", "sailfish"),
    txIn = TRUE,
    txOut = FALSE,
    tx2gene,
    ignoreTxVersion = TRUE
) {
    assert_all_are_dirs(sampleDirs)
    assert_are_identical(names(sampleDirs), basename(sampleDirs))
    type <- match.arg(type)
    assert_is_a_bool(txIn)
    assert_is_a_bool(txOut)
    assert_is_all_of(tx2gene, "tx2gene")
    tx2gene <- as.data.frame(tx2gene)

    # Locate the counts files --------------------------------------------------
    subdirs <- file.path(sampleDirs, type)
    assert_all_are_dirs(subdirs)

    if (type %in% c("salmon", "sailfish")) {
        # Use `quant.sf` file for salmon or sailfish.
        files <- file.path(subdirs, "quant.sf")
    } else if (type == "kallisto") {
        # Use `abundance.h5` file (HDF5) for kallisto.
        files <- file.path(subdirs, "abundance.h5")
    }
    assert_all_are_existing_files(files)
    names(files) <- names(sampleDirs)

    # Transcript versions ------------------------------------------------------
    # Ensure transcript IDs are stripped from tx2gene, if desired.
    if (isTRUE(ignoreTxVersion)) {
        tx2gene[["transcriptID"]] <-
            stripTranscriptVersions(tx2gene[["transcriptID"]])
        rownames(tx2gene) <- tx2gene[["transcriptID"]]
    }

    # Import counts using tximport ---------------------------------------------
    # Note that this step can take a long time when processing a lot of samples,
    # and is recommended to be run on an HPC cluster, rather than locally.
    message(paste(
        "Reading", type, "counts using tximport",
        packageVersion("tximport")
    ))
    txi <- tximport(
        files = files,
        type = type,
        txIn = txIn,
        txOut = txOut,
        countsFromAbundance = "lengthScaledTPM",
        tx2gene = tx2gene,
        ignoreTxVersion = ignoreTxVersion,
        importer = read_tsv
    )
    # Run assert checks on the matrices.
    invisible(lapply(
        X = txi[c(
            "abundance",
            "counts",
            "length"
        )],
        FUN = function(x) {
            assert_is_matrix(x)
            assert_are_identical(
                x = colnames(x),
                y = names(sampleDirs)
            )
        }
    ))

    # Ensure versions are stripped ---------------------------------------------
    # Counts from GENCODE annotations can require additional sanitization.
    if (isTRUE(ignoreTxVersion)) {
        genes <- rownames(txi[["abundance"]])
        genes <- stripTranscriptVersions(genes)
        assert_has_no_duplicates(genes)

        # Update the matrices.
        rownames(txi[["abundance"]]) <- genes
        rownames(txi[["counts"]]) <- genes
        rownames(txi[["length"]]) <- genes

        # Update the inferential replicates.
        infReps <- txi[["infReps"]]
        if (
            is.list(infReps) &&
            has_length(infReps)
        ) {
            infReps <- lapply(infReps, function(x) {
                rownames(x) <- genes
                x
            })
            txi[["infReps"]] <- infReps
        }
    }

    # Return -------------------------------------------------------------------
    .assertIsTximport(txi)
    txi
}



.assertIsTximport <- function(txi) {
    assert_is_list(txi)
    assert_are_intersecting_sets(
        x = c(
            "abundance",
            "counts",
            "infReps",  # v1.9+
            "length",
            "countsFromAbundance"
        ),
        y = names(txi)
    )

    abundance <- txi[["abundance"]]
    counts <- txi[["counts"]]
    infReps <- txi[["infReps"]]
    length <- txi[["length"]]
    countsFromAbundance <- txi[["countsFromAbundance"]]

    assert_are_identical(dimnames(abundance), dimnames(counts))
    assert_are_identical(dimnames(abundance), dimnames(length))

    # Inferential replicates added in v1.9.
    if (
        is.list(infReps) &&
        has_length(infReps)
    ) {
        assert_are_identical(names(infReps), colnames(abundance))
        assert_are_identical(rownames(infReps[[1L]]), rownames(abundance))
    }

    assert_is_a_string(countsFromAbundance)
    assert_is_non_empty(countsFromAbundance)
}



.regenerateTximport <- function(object) {
    validObject(object)
    assert_is_all_of(object, "RangedSummarizedExperiment")

    samples <- colnames(object)
    genes <- rownames(object)

    abundance <- assays(object)[["tpm"]]
    counts <- assays(object)[["counts"]]
    length <- assays(object)[["length"]]

    infReps <- metadata(object)[["infReps"]]
    if (is.list(infReps)) {
        infReps <- infReps[samples]
        infReps <- lapply(infReps, function(x) {
            x[genes, , drop = FALSE]
        })
    }

    countsFromAbundance <- metadata(object)[["countsFromAbundance"]]

    txi <- list(
        abundance = abundance,
        counts = counts,
        infReps = infReps,
        length = length,
        countsFromAbundance = countsFromAbundance
    )
    txi <- Filter(Negate(is.null), txi)

    .assertIsTximport(txi)
    txi
}
