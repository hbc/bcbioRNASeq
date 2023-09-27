#' @name plotPcaCovariates
#' @author Lorena Pantano, Michael Steinbaugh, Rory Kirchner
#' @inherit AcidGenerics::plotPcaCovariates
#' @note Requires the DEGreport package to be installed.
#' @note Updated 2022-10-24.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments, passed to `DEGreport::degCovariates()`.
#'
#' @param fdr `numeric(1)`.
#' Cutoff to determine the minimum false discovery rate (FDR) to consider
#' significant correlations between principal components (PCs) and covariates.
#'
#' @param metrics `boolean`. Include sample summary metrics as covariates.
#' Defaults to include all metrics columns (`TRUE`), but desired columns can
#' be specified here as a character vector.
#'
#' @seealso
#' - `DEGreport::degCovariates()`.
#' - `DESeq2::rlog()`.
#' - `DESeq2::varianceStabilizingTransformation()`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' if (requireNamespace("DEGreport", quietly = TRUE)) {
#'     plotPcaCovariates(bcb)
#' }
NULL



## Updated 2022-05-25.
`plotPcaCovariates,bcbioRNASeq` <- # nolint
    function(object,
             metrics = TRUE,
             normalized,
             fdr = 0.1) {
        assert(
            requireNamespaces("DEGreport"),
            validObject(object),
            isAny(metrics, c("character", "logical")),
            isNumber(fdr)
        )
        normalized <- match.arg(normalized)
        ## Get the normalized counts.
        counts <- counts(object, normalized = normalized)
        data <- metrics(object)
        ## Get factor (metadata) columns.
        keep <- which(bapply(data, is.factor))
        factors <- data[, keep, drop = FALSE]
        ## Drop columns that are all zeroes (not useful to plot). Sometimes we
        ## are missing values for some samples but not others; plotPcaCovariates
        ## was failing in those cases when checking if a column was a numeric.
        ## Here we ignore the NAs for numeric column checking. Adding the
        ## `na.rm` here fixes the issue.
        keep <- which(bapply(data, is.numeric))
        numerics <- data[, keep, drop = FALSE]
        keep <- which(colSums(as.data.frame(numerics), na.rm = TRUE) > 0L)
        numerics <- numerics[, keep, drop = FALSE]
        metadata <- cbind(factors, numerics)
        ## Select the metrics to use for plot.
        if (isTRUE(metrics)) {
            ## Sort columns alphabetically.
            col <- sort(colnames(metadata))
        } else if (identical(metrics, FALSE)) {
            ## Use the defined interesting groups.
            col <- interestingGroups(object)
        } else if (is.character(metrics)) {
            col <- metrics
        }
        ## Check for minimum number of metrics.
        if (length(col) < 2L) {
            abort(sprintf(
                "{.fun %s} requires >= 2 metrics.", "plotPcaCovariates"
            ))
        }
        assert(isSubset(col, colnames(metadata)))
        metadata <- metadata[, col]
        suppressWarnings({
            DEGreport::degCovariates(
                counts = counts,
                metadata = metadata,
                fdr = fdr
            )
        })
    }

formals(`plotPcaCovariates,bcbioRNASeq`)[["normalized"]] <- # nolint
    .normalized



#' @rdname plotPcaCovariates
#' @export
setMethod(
    f = "plotPcaCovariates",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotPcaCovariates,bcbioRNASeq`
)
