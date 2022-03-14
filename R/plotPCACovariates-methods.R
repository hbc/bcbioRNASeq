#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh, Rory Kirchner
#' @inherit AcidGenerics::plotPCACovariates
#' @note Requires the DEGreport package to be installed.
#' @note Updated 2022-03-07.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments, passed to `DEGreport::degCovariates()`.
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
#'     plotPCACovariates(bcb)
#' }
NULL



## Updated 2021-09-10.
`plotPCACovariates,bcbioRNASeq` <- # nolint
    function(object,
             metrics = TRUE,
             normalized,
             ...) {
        alertWarning(sprintf(
            "Disabled until bug is fixed in {.pkg %s}.", "DEGreport"
        ))
        return(invisible())
        validObject(object)
        assert(
            requireNamespace("DEGreport", quietly = TRUE),
            isAny(metrics, c("character", "logical"))
        )
        normalized <- match.arg(normalized)
        ## Get the normalized counts.
        counts <- counts(object, normalized = normalized)
        data <- metrics(object)
        ## Get factor (metadata) columns.
        keep <- which(bapply(data, is.factor))
        factors <- data[, keep, drop = FALSE]
        ## Drop columns that are all zeroes (not useful to plot). Sometimes we
        ## are missing values for some samples but not others; plotPCACovariates
        ## was failing in those cases when checking if a column was a numeric.
        ## Here we ignore the NAs for numeric column checking. Adding the
        ## `na.rm` here fixes the issue.
        keep <- which(bapply(data, is.numeric))
        numerics <- data[, keep, drop = FALSE]
        keep <- which(colSums(numerics, na.rm = TRUE) > 0L)
        numerics <- numerics[, keep, drop = FALSE]
        metadata <- cbind(factors, numerics)
        rownames(metadata) <- data[["sampleId"]]
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
                "{.fun %s} requires >= 2 metrics.", "plotPCACovariates"
            ))
        }
        assert(isSubset(col, colnames(metadata)))
        metadata <- metadata[, col]
        ## Handle warnings in DEGreport more gracefully.
        withCallingHandlers(
            expr = DEGreport::degCovariates(
                counts = counts,
                metadata = metadata,
                ...
            ),
            warning = function(w) {
                if (isTRUE(grepl(
                    pattern = "joining character vector and factor",
                    x = as.character(w)
                ))) {
                    invokeRestart("muffleWarning")
                } else if (isTRUE(grepl(
                    pattern = "Unquoting language objects",
                    x = as.character(w)
                ))) {
                    invokeRestart("muffleWarning")
                } else {
                    w
                }
            }
        )
    }

formals(`plotPCACovariates,bcbioRNASeq`)[["normalized"]] <- .normalized



#' @rdname plotPCACovariates
#' @export
setMethod(
    f = "plotPCACovariates",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotPCACovariates,bcbioRNASeq`
)
