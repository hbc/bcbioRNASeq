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
#' @note Updated 2021-09-10.
#'
#' @param sampleDirs `character`.
#'   Sample directories to import.
#' @param type `character(1)`.
#'   Expression caller to use.
#'
#' @seealso
#' - `tximport::tximport()`.
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
    basename <- ifelse(
        test = identical(type, "kallisto"),
        yes = "abundance.h5",
        ## salmon, sailfish.
        no = "quant.sf"
    )
    files <- file.path(subdirs, basename)
    assert(all(isFile(files)))
    names(files) <- names(sampleDirs)
    ## tx2gene -----------------------------------------------------------------
    tx2gene <- as.data.frame(tx2gene)
    if (!identical(
        x = c("txId", "txName"),
        y = colnames(tx2gene)
    )) {
        colnames(tx2gene) <- c("txId", "txName")
    }
    ## Ensure transcript IDs are stripped from tx2gene, if desired.
    if (isTRUE(ignoreTxVersion)) {
        tx2gene[["txId"]] <- stripTranscriptVersions(tx2gene[["txId"]])
        rownames(tx2gene) <- tx2gene[["txId"]]
    }
    ## Import counts using tximport --------------------------------------------
    ## Note that this step can take a long time when processing a lot of
    ## samples, and is recommended to be run on an HPC cluster, rather than
    ## locally.
    alert(sprintf(
        fmt = paste(
            "Importing {.pkg %s} transcript-level counts from {.file %s}",
            "files using {.pkg tximport} %s."
        ),
        type, basename(files[[1L]]), packageVersion("tximport")
    ))
    dl(c(
        "countsFromAbundance" = countsFromAbundance,
        "txOut" = txOut
    ))
    ## We're using a `do.call()` approach here so we can apply version-specific
    ## tximport fixes, if necessary.
    args <- list(
        "files" = files,
        "type" = type,
        "txIn" = TRUE,
        "txOut" = txOut,
        "countsFromAbundance" = countsFromAbundance,
        "tx2gene" = tx2gene,
        "ignoreTxVersion" = ignoreTxVersion,
        "geneIdCol" = "geneId",
        "txIdCol" = "txId"
    )
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



#' Detect if object is tximport list return
#'
#' @note Updated 2019-08-20.
#' @noRd
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
    assert(isString(countsFromAbundance))
    TRUE
}
