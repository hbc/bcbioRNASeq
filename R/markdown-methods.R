#' @name markdown
#' @importFrom basejump markdown
#' @inherit basejump::markdown
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams general
#'
#' @examples
#' data(deseq_small)
#' x <- DESeqResultsTables(deseq_small)
#' markdown(x)
NULL



# bcbioRNASeq ==================================================================
.markdown.bcbioRNASeq <-  # nolint
    function(object) {
        rse <- as(object, "RangedSummarizedExperiment")
        sampleData(rse) <- sampleData(object, clean = TRUE)
        markdown(rse)
    }



#' @describeIn markdown Sample metadata table.
#' @export
setMethod(
    f = "markdown",
    signature = signature("bcbioRNASeq"),
    definition = .markdown.bcbioRNASeq
)



# DESeqResultsTables ===========================================================
.markdown.DESeqResultsTables <-  # nolint
    function(
        object,
        headerLevel = 2L
    ) {
        stopifnot(is(object, "DESeqResultsTables"))
        validObject(object)
        assertIsHeaderLevel(headerLevel)
        
        metadata <- slot(object, "metadata")

        # Include a contrast header, which is useful for looping.
        contrast <- contrastName(object)
        markdownHeader(contrast, level = headerLevel, asis = TRUE)
        headerLevel <- headerLevel + 1L

        # File paths -----------------------------------------------------------
        files <- metadata[["files"]]

        if (length(files) > 0L) {
            # Get Dropbox URLs, if necessary.
            dropbox <- metadata[["dropbox"]]
            if (isTRUE(dropbox)) {
                # nocov start
                # Using local Dropbox token for code coverage here.
                files <- vapply(
                    X = files,
                    FUN = function(x) {
                        x[["url"]]
                    },
                    FUN.VALUE = "character"
                )
                # Remove trailing query string from URL.
                basenames <- gsub("\\?.*$", "", basename(files))
                # nocov end
            } else {
                basenames <- basename(files)
            }
            names(basenames) <- names(files)

            markdownHeader("File downloads", level = headerLevel, asis = TRUE)
            markdownList(c(
                paste0(
                    "[`", basenames[["all"]], "`]",
                    "(", files[["all"]], "): ",
                    "All genes, sorted by Ensembl identifier."
                ),
                paste0(
                    "[`", basenames[["deg"]], "`]",
                    "(", files[["deg"]], "): ",
                    "Genes that pass the alpha (FDR) and",
                    "log2 fold change (LFC) cutoffs."
                ),
                paste0(
                    "[`", basenames[["degUp"]], "`]",
                    "(", files[["degUp"]], "): ",
                    "Upregulated DEG; positive fold change."
                ),
                paste0(
                    "[`", basenames[["degDown"]], "`]",
                    "(", files[["degDown"]], "): ",
                    "Downregulated DEG; negative fold change."
                )
            ), asis = TRUE)
        } else {
            message("Object doesn't contain saved file paths.")
        }

        # Top tables -----------------------------------------------------------
        markdownHeader("Top tables", level = headerLevel, asis = TRUE)
        topTables(object)
    }



#' @describeIn markdown File paths (if exported) and top tables.
#' @export
setMethod(
    f = "markdown",
    signature = signature("DESeqResultsTables"),
    definition = .markdown.DESeqResultsTables
)



# DESeqAnalysis ================================================================
.markdown.DESeqAnalysis <-  # nolint
    function(object) {
        show(markdownHeader("Contrast names"))
        show(markdownList(.contrastNames(object)))
    }



#' @describeIn markdown List of contrast names.
#' @export
setMethod(
    f = "markdown",
    signature = signature("DESeqAnalysis"),
    definition = .markdown.DESeqAnalysis
)
