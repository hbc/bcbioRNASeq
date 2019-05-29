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
    type <- match.arg(type)
    assert_is_a_bool(txIn)
    assert_is_a_bool(txOut)
    tx2gene <- as.data.frame(tx2gene)
    assertIsTx2gene(tx2gene)

    # Check for count output format, by using the first sample directory.
    subdirs <- list.dirs(
        path = sampleDirs[[1L]],
        full.names = TRUE,
        recursive = FALSE
    )
    assert_are_intersecting_sets(basename(subdirs), validCallers)

    # Locate the counts files.
    subdirs <- file.path(sampleDirs, type)
    stopifnot(all(dir.exists(subdirs)))
    if (type %in% c("salmon", "sailfish")) {
        basename <- "quant.sf"
    } else if (type == "kallisto") {
        basename <- "abundance.h5"
    }
    files <- file.path(subdirs, basename)
    stopifnot(all(file.exists(files)))
    names(files) <- names(sampleDirs)

    # Begin loading of selected counts.
    message(paste(
        "Reading", type, "counts using tximport",
        packageVersion("tximport")
    ))

    # Ensure transcript IDs are stripped from tx2gene, if desired.
    if (isTRUE(ignoreTxVersion)) {
        tx2gene[["transcriptID"]] <-
            stripTranscriptVersions(tx2gene[["transcriptID"]])
        rownames(tx2gene) <- tx2gene[["transcriptID"]]
    }

    tximport(
        files = files,
        type = type,
        txIn = txIn,
        txOut = txOut,
        countsFromAbundance = "lengthScaledTPM",
        tx2gene = tx2gene,
        ignoreTxVersion = ignoreTxVersion,
        importer = read_tsv
    )
}



.regenerateTximportList <- function(object) {
    assert_is_any_of(
        x = object,
        classes = c("bcbioRNASeq", "RangedSummarizedExperiment")
    )
    assert_is_subset(c("tpm", "counts", "length"), names(assays(object)))
    assert_is_subset("countsFromAbundance", names(metadata(object)))
    list(
        abundance = assays(object)[["tpm"]],
        counts = assays(object)[["counts"]],
        length = assays(object)[["length"]],
        countsFromAbundance = metadata(object)[["countsFromAbundance"]]
    )
}
