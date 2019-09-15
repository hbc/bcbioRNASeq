## FIXME Gender markers plot needs a title by default.



#' @name plotGenderMarkers
#' @author Michael Steinbaugh
#' @inherit acidplots::plotGenderMarkers
#' @note Updated 2019-09-15.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in acidplots.
#'   See [acidplots::plotGenderMarkers()] for details.
#'
#' @examples
#' data(bcb)
#' rownames(bcb)[seq_len(2L)] <- c("ENSMUSG00000086503", "ENSMUSG00000069045")
#' plotGenderMarkers(bcb)
NULL



#' @rdname plotGenderMarkers
#' @name plotGenderMarkers
#' @importFrom bioverbs plotGenderMarkers
#' @importMethodsFrom acidplots plotGenderMarkers
#' @usage plotGenderMarkers(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotGenderMarkers,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        do.call(
            what = plotGenderMarkers,
            args = .dynamicTrans(
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
