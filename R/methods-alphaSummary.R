#' Print Summary Statistics of Alpha Level Cutoffs
#'
#' @note [bcbioRNASeq] does not support contrast definitions, since the
#'   object contains an internal [DESeqDataSet] with an empty design formula.
#'
#' @rdname alphaSummary
#' @name alphaSummary
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
#' @inheritParams DESeq2::results
#'
#' @param alpha Numeric vector of desired alpha cutoffs.
#' @param caption *Optional.* Character vector to add as caption to the table.
#' @param ... *Optional.* Passthrough arguments to [DESeq2::results()]. Use
#'   either `contrast` or `name` arguments to define the desired contrast.
#'
#' @return [kable].
#'
#' @seealso [DESeq2::results()].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioRNASeq"))
#' load(system.file("extdata/dds.rda", package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' alphaSummary(bcb)
#'
#' # DESeqDataSet
#' alphaSummary(dds)
#' alphaSummary(dds, contrast = c("group", "ko", "ctrl"))
#' alphaSummary(dds, name = "group_ko_vs_ctrl")
NULL



# Constructors =================================================================
#' @importFrom DESeq2 results
#' @importFrom dplyr bind_cols
#' @importFrom knitr kable
#' @importFrom magrittr set_colnames set_rownames
#' @importFrom utils capture.output
.alphaSummary <- function(
    dds,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    caption = NULL,
    ...) {
    assert_is_all_of(dds, "DESeqDataSet")
    assert_is_numeric(alpha)
    assert_is_a_string_or_null(caption)
    dots <- list(...)

    # Generate an automatic caption
    if (is.null(caption)) {
        if (!is.null(dots[["contrast"]])) {
            caption <- dots[["contrast"]] %>%
                paste(collapse = " ")
        } else if (!is.null(dots[["name"]])) {
            caption <- dots[["name"]]
        }
    }

    dflist <- lapply(seq_along(alpha), function(a) {
        output <- capture.output(summary(results(dds, ..., alpha = alpha[a])))
        # Subset the lines of interest from summary
        output <- output[4L:8L]
        # Extract the values after the colon in summary
        values <- vapply(
            X = output,
            FUN = function(a) {
                gsub("^.+\\:\\s(.+)\\s$", "\\1", a)
            },
            FUN.VALUE = "character")
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



#' @importFrom BiocGenerics design
#' @importFrom stats formula
.alphaSummary.bcbioRNASeq <- function(  # nolint
    object,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    caption = NULL,
    ...) {
    dds <- bcbio(object, "DESeqDataSet")
    # Warn if empty design formula detected
    if (design(dds) == formula(~1)) {  # nolint
        warn("Internal DESeqDataSet has an empty design formula")
    }
    .alphaSummary(
        dds = dds,
        alpha = alpha,
        caption = caption,
        ...)
}



.alphaSummary.DESeqDataSet <- function(  # nolint
    object,
    alpha = c(0.1, 0.05, 0.01, 1e-3, 1e-6),
    caption = NULL,
    ...) {
    .alphaSummary(
        dds = object,
        alpha = alpha,
        caption = caption,
        ...)
}



# Methods ======================================================================
#' @rdname alphaSummary
#' @export
setMethod(
    "alphaSummary",
    signature("bcbioRNASeq"),
    .alphaSummary.bcbioRNASeq)



#' @rdname alphaSummary
#' @export
setMethod(
    "alphaSummary",
    signature("DESeqDataSet"),
    .alphaSummary.DESeqDataSet)
