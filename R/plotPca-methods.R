#' Principal component analysis plot
#'
#' @name plotPca
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#' See `AcidPlots::plotPca()` for details.
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotPca(bcb, label = FALSE)
#' plotPca(bcb, label = TRUE)
NULL



## Updated 2020-09-15.
`plotPca,bcbioRNASeq` <- # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotPca,
            args = list(
                object = .normalizedSE(
                    object = object,
                    normalized = match.arg(normalized)
                ),
                ...
            )
        )
    }

formals(`plotPca,bcbioRNASeq`)[["normalized"]] <- # nolint
    .dt



#' @rdname plotPca
#' @export
setMethod(
    f = "plotPca",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotPca,bcbioRNASeq`
)
