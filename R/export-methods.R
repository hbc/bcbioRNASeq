#' @name export
#' @importFrom basejump export
#' @inherit basejump::export
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams bcbioBase::copyToDropbox
#' @inheritParams general
#' @param dir `string`. Local directory path.
#' @param dropboxDir `string` or `NULL`. Dropbox directory path where to archive
#'   the results tables for permanent storage (e.g. Stem Cell Commons). When
#'   this option is enabled, unique links per file are generated internally with
#'   the rdrop2 package. Note that local files are written to [base::tempdir()]
#'   and the `dir` argument is ignored, if this is enabled.
#' @param rdsToken `string` or `NULL`. RDS file token to use for Dropbox
#'   authentication. If set `NULL` and `dropboxDir` is defined, then an
#'   interactive prompt will appear requesting authorization.
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




.export.DESeqResultsTables <-  # nolint
    function(
        x,
        dir = ".",
        compress = TRUE,
        dropboxDir = NULL,
        rdsToken = NULL
    ) {
        validObject(x)
        # Write local files to tempdir if Dropbox mode is enabled.
        if (is_a_string(dropboxDir)) {
            dir <- tempdir()  # nocov
        } else {
            dir <- initializeDirectory(dir)
        }
        assert_is_a_bool(compress)
        assertIsAStringOrNULL(dropboxDir)

        # Extract the results tables from the object.
        tables <- coerceS4ToList(x)[c("all", "deg", "degUp", "degDown")]

        # Local files (required) -----------------------------------------------
        stem <- snake(contrastName(x))
        format <- "csv"
        if (isTRUE(compress)) {
            format <- paste0(format, ".gz")
        }
        ext <- paste0(".", format)
        localFiles <- file.path(
            dir,
            paste0(stem, "_", snake(names(tables)), ext)
        )
        names(localFiles) <- names(tables)

        # Write the results tables to local directory
        if (length(dropboxDir)) {
            # nocov start : use local Dropbox token
            message(paste0(
                "Writing ",
                toString(basename(localFiles)),
                " to Dropbox (", dropboxDir, ")."
            ))
            # nocov end
        } else {
            message(paste0(
                "Writing ", toString(basename(localFiles)),
                " to ", dir, "."
            ))
        }

        invisible(mapply(
            x = tables,
            file = localFiles,
            FUN = function(x, file) {
                export(x = x, file = file)
            },
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        ))

        # Check that all writes were successful.
        assert_all_are_existing_files(localFiles)

        # Assign the local file paths to the object.
        # Slot these only when we're not saving to Dropbox.
        if (is.null(dropboxDir)) {
            slot(x, "localFiles") <- localFiles
        }

        # Copy to Dropbox (optional) -------------------------------------------
        if (is.character(dropboxDir)) {
            # nocov start
            # Using a local token to cover this code.
            dropboxFiles <- copyToDropbox(
                files = localFiles,
                dir = dropboxDir,
                rdsToken = rdsToken
            )
            assert_is_list(dropboxFiles)
            slot(x, "dropboxFiles") <- dropboxFiles
            # nocov end
        }

        x
    }



#' @rdname export
#' @export
setMethod(
    f = "export",
    signature = signature("DESeqResultsTables"),
    definition = .export.DESeqResultsTables
)
