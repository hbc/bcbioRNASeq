# FIXME Rename the `use` argument, it's too vague
# FIXME Function is returning this error:
# Error in colnames(COVAR_data)[sapply(COVAR_data, is.factor)] :
# invalid subscript type 'list'



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
setMethod("plot_pca_covariates", "bcbioRNADataSet", function(
    object, transform = "rlog", use = NULL, ...) {
    metrics <- metrics(object)
    if (is.null(metrics)) return(NULL)

    # Check for valid `transform` argument
    transform_args <- c("rlog", "vst")
    if (!transform %in% transform_args) {
        stop(paste("Valid transforms:", toString(transform_args)))
    }

    # Get the columns of interest
    if (is.null(use)) {
        # FIXME We want to use `interesting_groups` by default right?
        use <- colnames(metrics)
    } else {
        use <- intersect(use, colnames(metrics))
    }
    if (length(use) == 0L) {
        stop("Not columns matched between use and metadata")
    }

    keep_metrics <- lapply(use, function(a) {
        if (length(unique(metrics[, a])) > 1L) a
    }) %>%
        unlist %>%
        .[!is.null(.)]

    metrics <- metrics %>%
        as.data.frame %>%
        set_rownames(.[["sample_name"]]) %>%
        .[, setdiff(keep_metrics, c("sample_name", "file_name")), drop = FALSE]

    # Pass internal [DESeqTransform] to [degCovariates()]
    if (!transform %in% c("rlog", "vst")) {
        stop("DESeqTransform must be rlog or vst")
    }

    res <- assays(object)[[transform]] %>%
        # Assay needed here to get the matrix from [DESeqTransform]
        assay %>%
        degCovariates(metadata = metrics, ...)
    # FIXME Need to clean up the plotting
    res %>%
        .[["plot"]] +
        theme(axis.text.x = element_text(angle = 60L, hjust = 1L))
    invisible(res)
})
