#' Plot Genes Detected
#'
#' @name plotGenesDetected
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @importFrom basejump plotGenesDetected
#' @inherit basejump::plotGenesDetected
#' @export
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotGenesDetected(bcb_small)
NULL



#' @rdname plotGenesDetected
#' @export
setMethod(
    f = "plotGenesDetected",
    signature = signature("bcbioRNASeq"),
    definition = getMethod("plotGenesDetected", "SummarizedExperiment")
)
