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
#' @param metrics Include sample summary metrics as covariates.
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
#' plotPCACovariates(bcb)
NULL



# Methods ====
#' @rdname plotPCACovariates
#' @export
setMethod("plotPCACovariates", "bcbioRNASeqANY", function(
    object,
    transform = "rlog",
    metrics = TRUE,
    ...) {
    # Check for valid transform argument
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
    if (!length(colnames(metadata))) {
        stop("Sample metadata is empty")
    }

    # Counts
    counts <- assays(object) %>%
        .[[transform]] %>%
        # Assay needed here to get the matrix from the slotted DESeqTransform
        assay()

    degCovariates(
        counts = counts,
        metadata = metadata,
        ...)
})
