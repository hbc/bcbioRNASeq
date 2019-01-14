#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNASeq] object.
#'
#' @name plotPCACovariates
#' @family Quality Control Functions
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams general
#' @param metrics `boolean`. Include sample summary metrics as covariates.
#'   Defaults to include all metrics columns (`TRUE`), but desired columns can
#'   be specified here as a character vector.
#' @param ... Additional arguments, passed to [DEGreport::degCovariates()].
#'
#' @seealso
#' - [DEGreport::degCovariates()].
#' - [DESeq2::rlog()].
#' - [DESeq2::varianceStabilizingTransformation()].
#'
#' @return `ggplot`.
#'
#' @examples
#' plotPCACovariates(
#'     object = bcb_small,
#'     metrics = c("exonicRate", "intronicRate")
#' )
NULL



#' @rdname plotPCACovariates
#' @export
setMethod(
    "plotPCACovariates",
    signature("bcbioRNASeq"),
    function(
        object,
        metrics = TRUE,
        normalized = c("vst", "rlog", "tmm", "tpm", "rle"),
        ...
    ) {
        validObject(object)
        assert_is_any_of(metrics, c("character", "logical"))
        normalized <- match.arg(normalized)

        # Get the normalized counts.
        counts <- counts(object, normalized = normalized)

        # Note that we're ungrouping the tibble here, otherwise the downstream
        # factor/numeric column selection won't work as expected.
        data <- ungroup(metrics(object))
        factors <- select_if(data, is.factor)
        # Drop columns that are all zeroes (not useful to plot).
        # @roryk Sometimes we are missing values for some samples but not
        # others; plotPCACovariates was failing in those cases when checking if
        # a column was a numeric. Here we ignore the NAs for numeric column
        # checking. Adding the `na.rm` here fixes the issue.
        numerics <- data %>%
            select_if(is.numeric) %>%
            .[, colSums(., na.rm = TRUE) > 0L, drop = FALSE]
        metadata <- cbind(factors, numerics)

        # Select the metrics to use for plot.
        if (identical(metrics, TRUE)) {
            # Sort columns alphabetically.
            col <- sort(colnames(metadata))
        } else if (identical(metrics, FALSE)) {
            # Use the defined interesting groups.
            col <- interestingGroups(object)
        } else if (is.character(metrics)) {
            col <- metrics
        }

        if (length(col) < 2L) {
            stop("plotPCACovariates() requires >= 2 metrics.")
        }
        assert_is_subset(col, colnames(metadata))
        metadata <- metadata[, col, drop = FALSE]

        # DEGreport v1.18 has issues joining character and factor columns.
        # - Warning: Column `compare`/`PC`
        # - Warning: Column `covar`/`term`
        suppressWarnings(
            degCovariates(counts = counts, metadata = metadata, ...)
        )
    }
)
