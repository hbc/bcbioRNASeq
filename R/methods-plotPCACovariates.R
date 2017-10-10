#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNASeq] object.
#'
#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param transform String specifying [DESeqTransform] slotted inside the
#'   [bcbioRNASeq] object:
#'   - `rlog` (**recommended**).
#'   - `vst`: variance stabilizing transformation.
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
#' data(bcb)
#' plotPCACovariates(bcb, metrics = TRUE)
#' plotPCACovariates(bcb, metrics = c("exonicRate", "intronicRate"))
#'
#' # Return `NULL` on example data, since only `sampleName` is interesting
#' plotPCACovariates(bcb, metrics = FALSE)
NULL



# Constructors ====
.plotPCACovariates <- function(
    object,
    transform = "rlog",
    metrics = TRUE,
    ...) {
    # Check for valid transform argument
    transformArgs <- c("rlog", "vst")
    if (!transform %in% transformArgs) {
        stop(paste("Valid transforms:", toString(transformArgs)))
    }

    metadata <- metrics(object)

    # Select the metrics to use for plot
    if (identical(metrics, TRUE)) {
        # Subset the metadata data.frame
        metadata <- metadata %>%
            # Select only the numeric columns
            select_if(is.numeric) %>%
            # Drop columns that are all zeroes (not useful to plot)
            .[, colSums(.) > 0]
        # Sort columns alphabetically
        metrics <- sort(colnames(metadata))
    } else if (identical(metrics, FALSE)) {
        # Use the defined interesting groups
        metrics <- interestingGroups(object)
    }

    # Now select the columns to use for plotting
    if (is.character(metrics) &
        all(metrics %in% colnames(metadata))) {
        metadata <- dplyr::select(metadata, metrics)
    } else {
        stop("Failed to select valid metrics for plot", call. = FALSE)
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



# Methods ====
#' @rdname plotPCACovariates
#' @export
setMethod(
    "plotPCACovariates",
    signature("bcbioRNASeqANY"),
    .plotPCACovariates)
