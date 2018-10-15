#' @inherit DESeqAnalysis-class
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @export
#'
#' @param data `DESeqDataSet`.
#' @param transform `DESeqTransform`.
#' @param results `list`. One or more unshrunken `DESeqResults`. Assign the
#'   [DESeq2::results()] return here.
#' @param lfcShrink `list`. One or more shrunken `DESeqResults`. Assign the
#'   [DESeq2::lfcShrink()] return here.
#'
#' @return `DESeqAnalysis`.
#'
#' @examples
#' data(bcb_small)
#' library(DESeq2)
#' dds <- as(bcb_small, "DESeqDataSet")
#' design(dds) <- ~ treatment
#' dds <- DESeq(dds)
#' class(dds)
#' vst <- varianceStabilizingTransformation(dds)
#' class(vst)
#' resultsNames(dds)
#' res <- results(dds, name = resultsNames(dds)[[2L]])
#' class(res)
#' x <- DESeqAnalysis(
#'     data = dds,
#'     transform = vst,
#'     results = list(res),
#'     lfcShrink = list(lfcShrink(dds = dds, coef = 2L))
#' )
#' print(x)
DESeqAnalysis <- function(
    data,
    transform,
    results,
    lfcShrink
) {
    new(
        Class = "DESeqAnalysis",
        data = data,
        transform = transform,
        results = results,
        lfcShrink = lfcShrink
    )
}



.contrastNames <- function(object) {
    vapply(
        X = object@results,
        FUN = contrastName,
        FUN.VALUE = character(1L)
    )
}
