#' @name plotPCACovariates
#' @author Lorena Pantano, Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotPCACovariates
#' @note Requires the DEGreport package to be installed.
#' @note Updated 2019-08-20.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
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
#' @examples
#' data(bcb)
#' plotPCACovariates(bcb)
NULL



#' @rdname plotPCACovariates
#' @name plotPCACovariates
#' @importFrom bioverbs plotPCACovariates
#' @usage plotPCACovariates(object, ...)
#' @export
NULL



## Updated 2019-08-20.
`plotPCACovariates,bcbioRNASeq` <-  # nolint
    function(
        object,
        metrics = TRUE,
        normalized,
        ...
    ) {
        validObject(object)
        assert(
            requireNamespace("DEGreport", quietly = TRUE),
            isAny(metrics, c("character", "logical"))
        )
        normalized <- match.arg(normalized)
        ## Get the normalized counts.
        counts <- counts(object, normalized = normalized)
        data <- metrics(object)
        factors <- data[, bapply(data, is.factor)]
        ## Drop columns that are all zeroes (not useful to plot). Sometimes we
        ## are missing values for some samples but not others; plotPCACovariates
        ## was failing in those cases when checking if a column was a numeric.
        ## Here we ignore the NAs for numeric column checking. Adding the
        ## `na.rm` here fixes the issue.
        numerics <- data[, which(bapply(data, is.numeric)), ]
        numerics <- numerics[, which(colSums(numerics, na.rm = TRUE) > 0L), ]
        metadata <- cbind(factors, numerics)
        rownames(metadata) <- data[["sampleID"]]
        ## Select the metrics to use for plot.
        if (identical(metrics, TRUE)) {
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
            stop("'plotPCACovariates()' requires >= 2 metrics.")
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

formals(`plotPCACovariates,bcbioRNASeq`)[["normalized"]] <- normalizedCounts



#' @rdname plotPCACovariates
#' @export
setMethod(
    f = "plotPCACovariates",
    signature = signature("bcbioRNASeq"),
    definition = `plotPCACovariates,bcbioRNASeq`
)
