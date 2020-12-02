#' @name plotPCA
#' @author Michael Steinbaugh
#' @importMethodsFrom AcidPlots plotPCA
#' @inherit AcidPlots::plotPCA
#' @note Updated 2020-09-15.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See [AcidPlots::plotPCA()] for details.
#'
#' @examples
#' data(bcb)
#' plotPCA(bcb, label = FALSE)
#' plotPCA(bcb, label = TRUE)
NULL



#' @rdname plotPCA
#' @name plotPCA
#' @importFrom BiocGenerics plotPCA
#' @usage plotPCA(object, ...)
#' @export
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
