# TODO Add dropbox mode for all supported methods.
# TODO Set `dropbox` as a boolean, instead of separate `dropboxDir` argument.
# TODO Need to add a bcbioRNASeq method that also exports the on-the-fly
# calculations (TMM, RLE).



#' @name export
#' @inherit basejump::export
#' @author Michael Steinbaugh
#'
#' @inheritParams bcbioBase::copyToDropbox
#' @inheritParams params
#'
#' @param counts `matrix`. Normalized counts. DESeq2 size-factor normalized
#'   counts or transcripts per million (TPM) are recommended.
#' @param dir `string`. Directory path.
#' @param dropbox `boolean`. Export results to [Dropbox][] instead of local
#'   storage. Recommended method by [HBC][] for permanent storage (e.g. [Stem
#'   Cell Commons][]). When this option is enabled, unique links per file are
#'   generated internally with [bcbioBase::copyToDropbox()], which relies on the
#'   [rdrop2][] package. Note that local files are written to [base::tempdir()]
#'   when this option is enabled.
#' @param rowData `DataFrame`. Row annotation data.
#' @param sampleNames Named `character`. Human readable sample names. Only
#'   applies when `counts` argument is defined. Names must correspond to
#'   `colnames` of `counts` (these should be valid in R; see
#'   [base::make.names()] for details). Values will be remapped onto the counts
#'   columns per sample in the exported file, and can contain non-alphanumeric
#'   characters, hyphens, spaces, or start with a number.
#'
#'   [Dropbox]: https://dropbox.com
#'   [HBC]: http://bioinformatics.sph.harvard.edu
#'   [rdrop2]: https://cran.r-project.org/package=rdrop2
#'   [Stem Cell Commons]: https://stemcellcommons.org
#' @param rdsToken `string` or `NULL`. RDS file token to use for Dropbox
#'   authentication. If set `NULL` and `dropbox = TRUE` then an interactive
#'   prompt will appear requesting authorization.
#'
#' @examples
#' data(deseq)
#'
#' ## DESeqResults ====
#' x <- as(deseq, "DESeqResults")
#' export(x, file = "example.csv")
#'
#' ## Clean up.
#' unlink("example.csv", recursive = TRUE)
#'
#' ## DESeqResultsTables ====
#' x <- DESeqResultsTables(deseq)
#' export(x, dir = "example")
#' list.files("example")
#'
#' ## Clean up.
#' unlink("example", recursive = TRUE)
NULL



#' @importFrom basejump export
#' @aliases NULL
#' @export
basejump::export



# Internal =====================================================================
.prepareDESeqResults <- function(
    object,
    rowData = NULL,
    counts = NULL,
    sampleNames = NULL
) {
    assert_is_all_of(object, "DESeqResults")
    assert_is_any_of(rowData, c("DataFrame", "NULL"))
    assert_is_any_of(counts, c("matrix", "NULL"))
    assert_is_any_of(sampleNames, c("character", "NULL"))

    # Coerce DESeqResults to DataFrame.
    data <- as(object, "DataFrame")

    # Row annotations.
    if (!is.null(rowData) && ncol(rowData) > 0L) {
        message("Joining annotations.")
        assert_is_all_of(rowData, "DataFrame")
        rowData <- sanitizeRowData(rowData)
        assert_are_identical(rownames(data), rownames(rowData))
        data <- cbind(data, rowData)
    }

    # Variance-stabilized counts (DESeqTransform).
    if (!is.null(counts) && ncol(counts) > 0L) {
        message("Joining counts.")
        assert_is_matrix(counts)
        assert_are_identical(rownames(data), rownames(counts))
        # Convert to human friendly sample names, if possible.
        if (
            is.character(sampleNames) &&
            has_length(sampleNames)
        ) {
            message("Mapping human-friendly sample names.")
            assert_has_names(sampleNames)
            assert_are_identical(names(sampleNames), colnames(counts))
            colnames(counts) <- as.character(sampleNames)
        }
        assert_are_disjoint_sets(colnames(data), colnames(counts))
        data <- cbind(data, counts)
    }

    data
}



# DESeqResults =================================================================
export.DESeqResults <-  # nolint
    function(
        x,
        file,
        format,
        rowData = NULL,
        counts = NULL,
        sampleNames = NULL
    ) {
        x <- .prepareDESeqResults(
            object = x,
            rowData = rowData,
            counts = counts,
            sampleNames = sampleNames
        )
        # Export using ANY method.
        export(x = x, file = file, format = format)
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature("DESeqResults"),
    definition = export.DESeqResults
)



# DESeqResultsTables ===========================================================
.prepareResultsTablesList <- function(object) {
    assert_that(is(object, "DESeqResultsTables"))
    validObject(object)

    deg <- slot(object, "deg")
    assert_is_list(deg)
    assert_are_identical(names(deg), c("up", "down"))

    # Up-regulated, down-regulated, and bidirectional DEGs.
    up <- deg[["up"]]
    down <- deg[["down"]]
    both <- c(up, down)

    results <- slot(object, "results")
    rowData <- slot(object, "rowRanges") %>%
        # Coerce to standard data.frame first, to collapse "X" ranges column.
        as.data.frame() %>%
        as("DataFrame")
    counts <- slot(object, "counts")
    sampleNames <- slot(object, "sampleNames")

    all <- .prepareDESeqResults(
        object = results,
        rowData = rowData,
        counts = counts,
        sampleNames = sampleNames
    )
    assert_that(is(all, "DataFrame"))

    list(
        all = all,
        up = all[up, , drop = FALSE],
        down = all[down, , drop = FALSE],
        both = all[both, , drop = FALSE]
    )
}



export.DESeqResultsTables <-  # nolint
    function(
        x,
        dir = ".",
        compress = FALSE,
        dropbox = FALSE,
        rdsToken = NULL
    ) {
        validObject(x)
        # Write local files to tempdir if Dropbox mode is enabled.
        if (isTRUE(dropbox)) {
            dir <- tempdir()  # nocov
        } else {
            dir <- initDir(dir)
        }
        assert_is_a_bool(compress)
        assert_is_a_bool(dropbox)

        # Prepare the subset tables.
        tables <- .prepareResultsTablesList(x)
        assert_is_list(tables)
        assert_are_identical(
            x = names(tables),
            y = c("all", "up", "down", "both")
        )

        # Local files (required) -----------------------------------------------
        stem <- snake(contrastName(x))
        format <- "csv"
        if (isTRUE(compress)) {
            format <- paste0(format, ".gz")
        }
        ext <- paste0(".", format)
        files <- file.path(
            dir,
            paste0(stem, "_", snake(names(tables)), ext)
        )
        names(files) <- names(tables)

        # Write the results tables to local directory.
        if (isTRUE(dropbox)) {
            # nocov start
            # Use local Dropbox token.
            message(paste0(
                "Writing ",
                toString(basename(files)),
                " to Dropbox (", dir, ")."
            ))
            # nocov end
        } else {
            message(paste0(
                "Writing ", toString(basename(files)), " to ", dir, "."
            ))
        }

        invisible(mapply(
            x = tables,
            file = files,
            FUN = function(x, file) {
                export(x = x, file = file)
            },
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        ))

        # Check that all writes were successful.
        assert_all_are_existing_files(files)

        # Copy to Dropbox (optional) -------------------------------------------
        if (isTRUE(dropbox)) {
            # nocov start
            # Using a local token to cover this code.
            files <- copyToDropbox(
                files = files,
                dir = dir,
                rdsToken = rdsToken
            )
            assert_is_list(files)
            slot(x, "metadata")[["dropbox"]] <- TRUE
            # nocov end
        }

        # Return ---------------------------------------------------------------
        # Assign the local file paths to the object.
        slot(x, "metadata")[["export"]] <- files

        x
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature("DESeqResultsTables"),
    definition = export.DESeqResultsTables
)
