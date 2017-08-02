#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plotPCACovariates
#'
#' @param transform String specifying [DESeqTransform] slotted inside the
#'   [bcbioRNADataSet]:
#'   - `rlog` (**recommended**).
#'   - `vst`: variance stabilizing transformation.
#' @param metrics Includesample summary metrics as covariates.
#' @param ... Additional arguments, passed to [degCovariates()].
#'
#' @seealso
#' - [DEGreport::degCovariates()].
#' - [DESeq2::rlog()].
#' - [DESeq2::varianceStabilizingTransformation()].
#'
#' @return [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plotPCACovariates(bcb)
setMethod("plotPCACovariates", "bcbioRNADataSet", function(
    object, transform = "rlog", metrics = TRUE, ...) {
    # Check for valid `transform` argument
    transformArgs <- c("rlog", "vst")
    if (!transform %in% transformArgs) {
        stop(paste("Valid transforms:", toString(transformArgs)))
    }

    # Metadata
    if (isTRUE(metrics)) {
        metadata <- metrics(object)
    } else {
        metadata <- .interestingColData(object)
    }
    metadata[["sampleName"]] <- NULL

    # Counts
    counts <- assays(object) %>%
        .[[transform]] %>%
        # Assay needed here to get the matrix from the slotted [DESeqTransform]
        assay

    degCovariates(counts, metadata)
})
