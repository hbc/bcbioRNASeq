#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNASeq] object.
#'
#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
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
.plotPCACovariates <- function(
    object,
    metrics = TRUE,
    normalized = "rlog",
    ...) {
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
    } else {
        abort("`metrics` must be a logical or character vector")
    }

    # Stop on 1 column
    if (length(col) == 1L) {
        abort(paste(
            "`degCovariates()` requires at least 2 metadata columns"
        ))
    }

    # Now select the columns to use for plotting
    if (all(col %in% colnames(metadata))) {
        metadata <- metadata[, col, drop = FALSE]
    } else {
        # FIXME Make this message more informative. Which metrics?
        abort("Failed to select valid `metrics` columns to plot")
    }

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
    .plotPCACovariates)
