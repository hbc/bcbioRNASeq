#' @name export
#' @importFrom basejump export
#' @inherit basejump::export
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams bcbioBase::copyToDropbox
#' @inheritParams general
#' @param dir `string`. Directory path.
#' @param dropbox `boolean`. Export results to [Dropbox][] instead of local
#'   storage. Recommended method by [HBC][] for permanent storage (e.g. [Stem
#'   Cell Commons][]). When this option is enabled, unique links per file are
#'   generated internally with [bcbioBase::copyToDropbox()], which relies on the
#'   [rdrop2][] package. Note that local files are written to [base::tempdir()]
#'   when this option is enabled.
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
#' data(deseq_small)
#' x <- DESeqResultsTables(deseq_small)
#' export(x, dir = "example")
#' list.files("example")
#'
#' ## Clean up.
#' unlink("example", recursive = TRUE)
NULL



# FIXME Add Dropbox support for `bcbioRNASeq`
# FIXME Also need to export the on-the-fly calculations (TMM, RLE).
# TODO Add dropbox mode for bcbioRNASeq, similar to DESeqResultsTables.
# Set `dropbox` as a boolean, instead of a separate dir call.



.resultsTablesList <- function(object) {
    stopifnot(is(object, "DESeqResultsTables"))
    validObject(object)

    results <- slot(object, "results")
    stopifnot(is(results, "DESeqResults"))

    deg <- slot(object, "deg")
    assert_is_list(deg)
    assert_are_identical(names(deg), c("up", "down"))

    up <- deg[["up"]]
    down <- deg[["down"]]
    both <- c(up, down)

    list(
        all = results,
        up = results[up, , drop = FALSE],
        down = results[down, , drop = FALSE],
        both = results[both, , drop = FALSE]
    )
}



# FIXME Need to update this method after we change the DESeqResultsTables structure.
# FIXME Ensure we're joining rowData and counts automatically here.
.export.DESeqResultsTables <-  # nolint
    function(
        x,
        dir = ".",
        compress = TRUE,
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

        tables <- .resultsTablesList(x)
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
                "Writing ", toString(basename(files)),
                " to ", dir, "."
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
                dir = dropboxDir,
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
    definition = .export.DESeqResultsTables
)
