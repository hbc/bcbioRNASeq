#' @name plotGenderMarkers
#' @author Michael Steinbaugh
#' @inherit AcidPlots::plotGenderMarkers
#' @note Updated 2019-09-16.
#'
#' @inheritParams plotCounts
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in AcidPlots.
#'   See [AcidPlots::plotGenderMarkers()] for details.
#'
#' @examples
#' data(bcb)
#' rownames(bcb)[seq_len(2L)] <- c("ENSMUSG00000086503", "ENSMUSG00000069045")
#' plotGenderMarkers(bcb)
NULL



#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @importFrom AcidGenerics plotGenderMarkers
#' @importMethodsFrom AcidPlots plotGenderMarkers
#' @usage plotGenderMarkers(object, ...)
#' @export
NULL



## Updated 2019-09-16.
`plotGenderMarkers,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotGenderMarkers,
            args = .normalizedPlotArgs(
                object = object,
                normalized = match.arg(normalized),
                ...
            )
        )
    }

formals(`plotGenderMarkers,bcbioRNASeq`)[["normalized"]] <- .normalized



#' @rdname plotGenderMarkers
#' @export
setMethod(
    f = "plotGenderMarkers",
    signature = signature("bcbioRNASeq"),
    definition = `plotGenderMarkers,bcbioRNASeq`
)
