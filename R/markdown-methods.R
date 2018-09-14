#' Markdown
#'
#' Markdownify an S4 object.
#'
#' @name markdown
#' @author Michael Steinbaugh
#' @importFrom basejump markdown
#' @export
#'
#' @return [writeLines()] output.
#'
#' @examples
#' object <- resultsTables(deseq_small)
#' markdown(object)
NULL



#' @rdname markdown
#' @export
setMethod(
    "markdown",
    signature("DESeqResultsTables"),
    function(
        object,
        headerLevel = 2L
    ) {
        assert_is_all_of(object, "DESeqResultsTables")
        validObject(object)
        assertIsImplicitInteger(headerLevel)

        # Include a contrast header, which is useful for looping.
        contrast <- contrastName(object)
        markdownHeader(contrast, level = headerLevel, asis = TRUE)
        headerLevel <- headerLevel + 1L

        # CSV file links -------------------------------------------------------
        paths <- NULL
        # Prioritze `dropboxFiles` over `localFiles`.
        if (length(object@dropboxFiles) > 0L) {
            # nocov start
            # Using local Dropbox token for code coverage here.
            paths <- vapply(
                X = object@dropboxFiles,
                FUN = function(x) {
                    x[["url"]]
                },
                FUN.VALUE = "character"
            )
            basenames <- gsub("\\?.*$", "", basename(paths))
            # nocov end
        } else if (length(object@localFiles) > 0L) {
            paths <- object@localFiles
        } else {
            message("Object doesn't contain saved file paths")
        }

        if (length(paths) > 0L) {
            basenames <- basename(paths)
            names(basenames) <- names(paths)
            markdownHeader("CSV file links", level = headerLevel, asis = TRUE)
            markdownList(c(
                paste0(
                    "[`", basenames[["all"]], "`]",
                    "(", paths[["all"]], "): ",
                    "All genes, sorted by Ensembl identifier."
                ),
                paste0(
                    "[`", basenames[["deg"]], "`]",
                    "(", paths[["deg"]], "): ",
                    "Genes that pass the alpha (FDR) and",
                    "log2 fold change (LFC) cutoffs."
                ),
                paste0(
                    "[`", basenames[["degUp"]], "`]",
                    "(", paths[["degUp"]], "): ",
                    "Upregulated DEG; positive fold change."
                ),
                paste0(
                    "[`", basenames[["degDown"]], "`]",
                    "(", paths[["degDown"]], "): ",
                    "Downregulated DEG; negative fold change."
                )
            ), asis = TRUE)
        }

        # Top tables -----------------------------------------------------------
        markdownHeader("Top tables", level = headerLevel, asis = TRUE)
        topTables(object)
    }

)
