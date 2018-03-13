#' Alpha Level Cutoff Summary Statistics
#'
#' @note `bcbioRNASeq` class does currently support contrast definitions, since
#'   the object contains an internal `DESeqDataSet` with an empty design
#'   formula.
#'
#' @name alphaSummary
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
#' @inheritParams DESeq2::results
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param caption *Optional.* Character vector to add as caption to the table.
#' @param ... *Optional.* Passthrough arguments to [DESeq2::results()]. Use
#'   either `contrast` or `name` arguments to define the desired contrast.
#'
#' @return `kable`.
#'
#' @seealso [DESeq2::results()].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq ====
#' alphaSummary(bcb)
#'
#' # DESeqDataSet ====
#' resultsNames(dds)
#' alphaSummary(dds)
#' alphaSummary(dds, contrast = c("day", "7", "0"))
#' alphaSummary(dds, name = "day_7_vs_0")
NULL



# Methods ======================================================================
#' @rdname alphaSummary
#' @importFrom DESeq2 results
#' @importFrom dplyr bind_cols
#' @importFrom knitr kable
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom utils capture.output
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

        # Abort on empty design formula
        if (design(object) == ~ 1) {  # nolint
            abort("Empty design formula detected")
        }

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



#' @rdname alphaSummary
#' @importFrom BiocGenerics design
#' @export
setMethod(
    "alphaSummary",
    signature("bcbioRNASeq"),
    function(
        object,
        alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
        caption = NULL,
        ...
    ) {
        validObject(object)
        dds <- assays(object)[["dds"]]
        assert_is_all_of(dds, "DESeqDataSet")
        alphaSummary(
            object = dds,
            alpha = alpha,
            caption = caption,
            ...
        )
    }
)
