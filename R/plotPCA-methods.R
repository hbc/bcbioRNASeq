#' Principal component analysis plot
#'
#' @name plotPCA
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#' See `AcidPlots::plotPCA()` for details.
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotPCA(bcb, label = FALSE)
#' plotPCA(bcb, label = TRUE)
NULL



## Updated 2020-09-15.
`plotPCA,bcbioRNASeq` <- # nolint
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
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotPCA,bcbioRNASeq`
)
