#' Import transcript-level abundances and counts
#'
#' Import RNA-seq counts using
#' [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport).
#'
#' Normalized counts are loaded as length-scaled transcripts per million.
#' Consult this [vignette](https://goo.gl/h6fm15) for more information.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @keywords internal
#' @noRd
#'
#' @note Ignoring transcript versions should work by default. There may be
#' an issue with genomes containing non-Ensembl transcript IDs, such as
#' C. elegans, although we need to double check.
#' @note Updated 2019-08-07.
#'
#' @inheritParams tximport::tximport
#' @param sampleDirs `character`.
#'   Sample directories to import.
#' @param type `character(1)`.
#'   Expression caller to use.
#'
#' @seealso [tximport::tximport()].
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
        identical(names(sampleDirs), makeNames(basename(sampleDirs))),
        isFlag(txOut),
        is(tx2gene, "Tx2Gene")
    )
    type <- match.arg(type)
    countsFromAbundance <- match.arg(
        arg = countsFromAbundance,
        choices = eval(formals(tximport)[["countsFromAbundance"]])
    )

    ## Locate the counts files -------------------------------------------------
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

    ## tx2gene -----------------------------------------------------------------
    tx2gene <- as.data.frame(tx2gene)
    if (isTRUE(ignoreTxVersion)) {
        ## Ensure transcript IDs are stripped from tx2gene.
        tx2gene[["transcriptID"]] <-
            stripTranscriptVersions(tx2gene[["transcriptID"]])
        rownames(tx2gene) <- tx2gene[["transcriptID"]]
    }

    ## Import counts using tximport --------------------------------------------
    ## Note that this step can take a long time when processing a lot of
    ## samples, and is recommended to be run on an HPC cluster, rather than
    ## locally.
    message(sprintf(
        "Reading %s transcript-level counts from %s files using tximport %s.",
        type, basename(files[[1L]]), packageVersion("tximport")
    ))
    if (countsFromAbundance != "no") {
        message(sprintf("Scaling transcripts using %s.", countsFromAbundance))
    }
    if (identical(txOut, FALSE)) {
        message("Returning transcript abundance at gene level.")
    }

    args <- list(
        files = files,
        type = type,
        txIn = TRUE,
        txOut = txOut,
        countsFromAbundance = countsFromAbundance,
        tx2gene = tx2gene,
        ignoreTxVersion = ignoreTxVersion
    )

    ## tximport version-specific settings.
    ## This ensures backward compatibility with R 3.4.
    if (packageVersion("tximport") < "1.10") {
        args[["importer"]] <- readr::read_tsv
    } else {
        ## Defaults to `read_tsv()` if installed, don't need to set.
        args[["importer"]] <- NULL
        args[["existenceOptional"]] <- FALSE
        args[["ignoreAfterBar"]] <- FALSE
        ## Ensure that we're mapping these IDs correctly.
        args[["geneIdCol"]] <- "geneID"
        args[["txIdCol"]] <- "transcriptID"
    }

    txi <- do.call(what = tximport, args = args)

    ## Assert checks before return.
    invisible(lapply(
        X = txi[c("abundance", "counts", "length")],
        FUN = function(x) {
            assert(
                is.matrix(x),
                identical(colnames(x), names(sampleDirs))
            )
        }
    ))
    assert(.isTximportReturn(txi))

    txi
}



## Detect if object is tximport list return.
## Updated 2019-07-23.
.isTximportReturn <- function(list) {
    assert(
        is.list(list),
        areIntersectingSets(
            x = c(
                "abundance",
                "counts",
                "infReps",  # v1.9+
                "length",
                "countsFromAbundance"
            ),
            y = names(list)
        )
    )

    abundance <- list[["abundance"]]
    counts <- list[["counts"]]
    infReps <- list[["infReps"]]
    length <- list[["length"]]
    countsFromAbundance <- list[["countsFromAbundance"]]

    assert(
        identical(dimnames(abundance), dimnames(counts)),
        identical(dimnames(abundance), dimnames(length))
    )

    ## Inferential replicates added in v1.9.
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
