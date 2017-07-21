#' Find Correlation Between Principal Components (PCs) and Covariates
#'
#' [DEGreport::degCovariates()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plot_pca_covariates
#' @author Lorena Pantano, Michael Steinbaugh
#'
#' @param transform String specifying [DESeqTransform] slotted inside the
#'   [bcbioRNADataSet]:
#'   - `rlog` (**recommended**).
#'   - `vst`: variance stabilizing transformation.
#' @param use *Optional*. Character vector of columns to use.
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
#' plot_pca_covariates(bcb)
setMethod("plot_pca_covariates", "bcbioRNADataSet", function(
    object, transform = "rlog", use = NULL, ...) {
    # Check for valid `transform` argument
    transform_args <- c("rlog", "vst")
    if (!transform %in% transform_args) {
        stop(paste("Valid transforms:", toString(transform_args)))
    }

    # Subset the metadata with the `use` column, if desired
    metadata <- .interesting_col_data(object)
    if (!is.null(use)) {
        metadata <- metadata[, use]
    }

    counts <- assays(object) %>%
        .[[transform]] %>%
        # Assay needed here to get the matrix from the slotted [DESeqTransform]
        assay

    degCovariates(counts, metadata)
})
