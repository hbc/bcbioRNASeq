#' Import transcript-level abundances and counts
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
#' @param sampleDirs `character`.
#'   Sample directories to import.
#' @param type `character(1)`.
#'   Expression caller to use.
#'
#' @seealso `tximport::tximport`.
#'
#' @return `list`.
.tximport <- function(
    sampleDirs,
    type = c("salmon", "kallisto", "sailfish"),
    txOut = FALSE,
    countsFromAbundance = "lengthScaledTPM",
    tx2gene,
    ignoreTxVersion = TRUE
) {
    assert(
        all(isDirectory(sampleDirs)),
        identical(names(sampleDirs), basename(sampleDirs)),
        isFlag(txOut),
        is(tx2gene, "Tx2Gene")
    )
    type <- match.arg(type)
    countsFromAbundance <- match.arg(
        arg = countsFromAbundance,
        choices = eval(formals(tximport)[["countsFromAbundance"]])
    )

    # Locate the counts files --------------------------------------------------
    subdirs <- file.path(sampleDirs, type)
    assert(all(isDirectory(subdirs)))
    if (type %in% c("salmon", "sailfish")) {
        basename <- "quant.sf"
    } else if (type == "kallisto") {
        basename <- "abundance.h5"
    }
    files <- file.path(subdirs, basename)
    assert(all(isFile(files)))
    names(files) <- names(sampleDirs)

    # tx2gene ------------------------------------------------------------------
    tx2gene <- as.data.frame(tx2gene)
    if (isTRUE(ignoreTxVersion)) {
        # Ensure transcript IDs are stripped from tx2gene.
        tx2gene[["transcriptID"]] <-
            stripTranscriptVersions(tx2gene[["transcriptID"]])
        rownames(tx2gene) <- tx2gene[["transcriptID"]]
    }

    # Import counts using tximport ---------------------------------------------
    # Note that this step can take a long time when processing a lot of samples,
    # and is recommended to be run on an HPC cluster, rather than locally.
    message(paste0(
        "Reading ", type, " counts using tximport ",
        packageVersion("tximport"), "."
    ))
    message(paste("Reading from", basename(files[[1L]]), "files."))
    if (countsFromAbundance != "no") {
        message(paste0("Scaling using ", countsFromAbundance, "."))
    }
    txi <- tximport(
        files = files,
        type = type,
        txIn = TRUE,
        txOut = txOut,
        countsFromAbundance = countsFromAbundance,
        tx2gene = tx2gene,
        ignoreTxVersion = ignoreTxVersion,
        importer = read_tsv
    )

    # Assert checks before return.
    invisible(lapply(
        X = txi[c(
            "abundance",
            "counts",
            "length"
        )],
        FUN = function(x) {
            assert(
                is.matrix(x),
                identical(colnames(x), names(sampleDirs))
            )
        }
    ))
    assert(.isTximport(txi))

    txi
}



.isTximport <- function(txi) {
    assert(
        is.list(txi),
        areIntersectingSets(
            x = c(
                "abundance",
                "counts",
                "infReps",  # v1.9+
                "length",
                "countsFromAbundance"
            ),
            y = names(txi)
        )
    )

    abundance <- txi[["abundance"]]
    counts <- txi[["counts"]]
    infReps <- txi[["infReps"]]
    length <- txi[["length"]]
    countsFromAbundance <- txi[["countsFromAbundance"]]

    assert(
        identical(dimnames(abundance), dimnames(counts)),
        identical(dimnames(abundance), dimnames(length))
    )

    # Inferential replicates added in v1.9.
    if (is.list(infReps) && hasLength(infReps)) {
        assert(
            identical(names(infReps), colnames(abundance)),
            identical(rownames(infReps[[1L]]), rownames(abundance))
        )
    }

    assert(
        isString(countsFromAbundance),
        isNonEmpty(countsFromAbundance)
    )

    TRUE
}
