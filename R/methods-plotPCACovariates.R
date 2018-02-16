#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNASeq] object.
#'
#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotPCA
#'
#' @param metrics Include sample summary metrics as covariates. Defaults to
#'   include all metrics columns (`TRUE`), but desired columns can be specified
#'   here as a character vector.
#' @param ... Additional arguments, passed to [DEGreport::degCovariates()].
#'
#' @seealso
#' - [DEGreport::degCovariates()].
#' - [DESeq2::rlog()].
#' - [DESeq2::varianceStabilizingTransformation()].
#'
#' @return [ggplot].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotPCACovariates(bcb, metrics = c("exonicRate", "intronicRate"))
NULL



# Constructors =================================================================
#' @importFrom DEGreport degCovariates
#' @importFrom dplyr select_if
.plotPCACovariates.bcbioRNASeq <- function(  # nolint
    object,
    metrics = TRUE,
    normalized = "rlog",
    ...) {
    assert_is_any_of(metrics, c("character", "logical"))
    assert_is_a_string(normalized)

    counts <- counts(object, normalized = normalized)

    metadata <- metrics(object)
    factors <- select_if(metadata, is.factor)
    numerics <- select_if(metadata, is.numeric) %>%
        # Drop columns that are all zeroes (not useful to plot)
        .[, colSums(.) > 0L]
    metadata <- cbind(factors, numerics)

    # Select the metrics to use for plot
    if (identical(metrics, TRUE)) {
        # Sort columns alphabetically
        col <- sort(colnames(metadata))
    } else if (identical(metrics, FALSE)) {
        # Use the defined interesting groups
        col <- interestingGroups(object)
    } else if (is.character(metrics)) {
        col <- metrics
    }

    # Stop on 1 metrics column
    assert_all_are_greater_than(length(col), 1L)
    assert_is_subset(col, colnames(metadata))
    metadata <- metadata[, col, drop = FALSE]

    degCovariates(
        counts = counts,
        metadata = metadata,
        ...)
}



# Methods ======================================================================
#' @rdname plotPCACovariates
#' @export
setMethod(
    "plotPCACovariates",
    signature("bcbioRNASeq"),
    .plotPCACovariates.bcbioRNASeq)
