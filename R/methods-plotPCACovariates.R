#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNASeq] object.
#'
#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param metrics Include sample summary metrics as covariates. Defaults to
#'   include all metrics columns (`TRUE`), but desired columns can be specified
#'   here as a character vector.
#' @param transform String specifying [DESeqTransform] slotted inside the
#'   [bcbioRNASeq] object:
#'   - `rlog` (**recommended**).
#'   - `vst`: variance stabilizing transformation.
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
    transform = "rlog",
    ...) {
    # Check for valid transform argument
    transformArgs <- c("rlog", "vst")
    if (!transform %in% transformArgs) {
        stop(paste("Valid transforms:", toString(transformArgs)))
    }

    metadata <- metrics(object)
    factors <- select_if(metadata, is.factor)
    numerics <- select_if(metadata, is.numeric) %>%
        # Drop columns that are all zeroes (not useful to plot)
        .[, colSums(.) > 0]
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
        stop("'metrics' must be 'TRUE/FALSE' or character vector",
             call. = FALSE)
    }

    # Stop on 1 column
    if (length(col) == 1) {
        stop(paste(
            "'degCovariates()' requires at least 2 metadata columns"
        ), call. = FALSE)
    }

    # Now select the columns to use for plotting
    if (all(col %in% colnames(metadata))) {
        metadata <- metadata[, col, drop = FALSE]
    } else {
        stop("Failed to select valid 'metrics' for plot", call. = FALSE)
    }

    # Counts
    counts <- assays(object) %>%
        .[[transform]] %>%
        # Assay needed to get matrix from the slotted DESeqTransform
        assay()

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
