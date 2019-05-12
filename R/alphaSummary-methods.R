#' Alpha Level Cutoff Summary Statistics
#'
#' Quickly generate a summary table of various alpha level cutoffs, for use in
#' an R Markdown report.
#'
#' @note `bcbioRNASeq` class does currently support contrast definitions, since
#'   the object contains an internal `DESeqDataSet` with an empty design
#'   formula.
#'
#' @name alphaSummary
#' @family Differential Expression Functions
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
#' @inheritParams DESeq2::results
#' @param alpha `numeric`. Multiple alpha cutoffs.
#' @param caption `string` or `NULL`. Table caption. If set `NULL`, will be
#'   generated automatically from the contrast used.
#' @param ... Passthrough arguments to [DESeq2::results()]. Use either
#'   `contrast` or `name` arguments to define the desired contrast.
#'
#' @return `kable`.
#'
#' @seealso [DESeq2::results()].
#'
#' @examples
#' # DESeqDataSet ====
#' design(dds_small)
#' resultsNames(dds_small)
#' alphaSummary(dds_small, contrast = c("treatment", "folic_acid", "control"))
#' alphaSummary(dds_small, name = "treatment_folic_acid_vs_control")
NULL



#' @rdname alphaSummary
#' @name alphaSummary
#' @importFrom bioverbs alphaSummary
#' @usage alphaSummary(object, ...)
#' @export
NULL



#' @rdname alphaSummary
#' @export
setMethod(
    "alphaSummary",
    signature("DESeqDataSet"),
    function(
        object,
        alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
        caption = NULL,
        ...
    ) {
        validObject(object)
        assert_is_numeric(alpha)
        assertIsAStringOrNULL(caption)
        dots <- list(...)

        # Generate an automatic caption
        if (is.null(caption)) {
            if (!is.null(dots[["contrast"]])) {
                caption <- paste(dots[["contrast"]], collapse = " ")
            } else if (!is.null(dots[["name"]])) {
                caption <- dots[["name"]]
            }
        }

        dflist <- lapply(alpha, function(x) {
            output <- capture.output(
                summary(results(object, ..., alpha = x))
            )
            # Subset the lines of interest from summary
            output <- output[4L:8L]
            # Extract the values after the colon in summary
            values <- vapply(
                X = output,
                FUN = function(x) {
                    gsub("^.+\\:\\s(.+)\\s$", "\\1", x)
                },
                FUN.VALUE = "character"
            )
            data.frame(alpha = values)
        })

        dflist %>%
            bind_cols() %>%
            set_colnames(alpha) %>%
            set_rownames(c(
                "LFC > 0 (up)",
                "LFC < 0 (down)",
                "outliers",
                "low counts",
                "cutoff"
            )) %>%
            kable(caption = caption)
    }
)
