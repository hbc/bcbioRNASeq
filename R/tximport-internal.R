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
    list <- tximport(
        files = files,
        type = type,
        txIn = txIn,
        txOut = txOut,
        countsFromAbundance = "lengthScaledTPM",
        tx2gene = tx2gene,
        ignoreTxVersion = ignoreTxVersion,
        importer = read_tsv
    )

    # Assert checks on tximport list -------------------------------------------
    assert_is_list(list)
    assert_are_identical(
        names(list),
        c("abundance", "counts", "length", "countsFromAbundance")
    )
    # Run assert checks on the matrices.
    invisible(lapply(
        X = list[c("abundance", "counts", "length")],
        FUN = function(x) {
            assert_is_matrix(x)
            assert_are_identical(
                x = colnames(x),
                y = names(sampleDirs)
            )
        }
    ))
    assert_is_non_empty(list[["countsFromAbundance"]])

    # Ensure versions are stripped ---------------------------------------------
    # Counts from GENCODE annotations can require additional sanitization.
    if (isTRUE(ignoreTxVersion)) {
        rownames <- rownames(list[["abundance"]])
        assert_are_identical(rownames, rownames(list[["counts"]]))
        assert_are_identical(rownames, rownames(list[["length"]]))
        rownames <- stripTranscriptVersions(rownames)
        assert_has_no_duplicates(rownames)
        rownames(list[["abundance"]]) <- rownames
        rownames(list[["counts"]]) <- rownames
        rownames(list[["length"]]) <- rownames
    }

    # Return -------------------------------------------------------------------
    list
}



.regenerateTximportList <- function(object) {
    assert_is_all_of(object, "RangedSummarizedExperiment")
    assert_is_subset(c("tpm", "counts", "length"), assayNames(object))
    assert_is_subset("countsFromAbundance", names(metadata(object)))
    list(
        abundance = assays(object)[["tpm"]],
        counts = assays(object)[["counts"]],
        length = assays(object)[["length"]],
        countsFromAbundance = metadata(object)[["countsFromAbundance"]]
    )
}
