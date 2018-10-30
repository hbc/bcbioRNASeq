#' @inherit DESeqAnalysis-class
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
#' library(DESeq2)
#' data(bcb)
#' 
#' dds <- as(bcb, "DESeqDataSet")
#' design(dds) <- ~ treatment
#' dds <- DESeq(dds)
#' class(dds)
#'
#' vst <- varianceStabilizingTransformation(dds)
#' class(vst)
#'
#' resultsNames(dds)
#' res <- results(dds, name = resultsNames(dds)[[2L]])
#' class(res)
#'
#' x <- DESeqAnalysis(
#'     data = dds,
#'     transform = vst,
#'     results = list(res),
#'     lfcShrink = list(lfcShrink(dds = dds, coef = 2L))
#' )
#' print(x)
DESeqAnalysis <-  # nolint
    function(
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
    assert_that(is(object, "DESeqAnalysis"))
    vapply(
        X = slot(object, "results"),
        FUN = contrastName,
        FUN.VALUE = character(1L)
    )
}
