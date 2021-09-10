#' @name plotPCA
#' @author Michael Steinbaugh
#' @importMethodsFrom AcidPlots plotPCA
#' @inherit AcidGenerics::plotPCA
#' @note Updated 2021-09-10.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See `AcidPlots::plotPCA()` for details.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotPCA(bcb, label = FALSE)
#' plotPCA(bcb, label = TRUE)
NULL



## Updated 2020-09-15.
`plotPCA,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotPCA,
            args = list(
                object = .normalizedSE(
                    object = object,
                    normalized = match.arg(normalized)
                ),
                ...
            )
        )
    }

formals(`plotPCA,bcbioRNASeq`)[["normalized"]] <- .dt



#' @rdname plotPCA
#' @export
setMethod(
    f = "plotPCA",
    signature = signature("bcbioRNASeq"),
    definition = `plotPCA,bcbioRNASeq`
)
