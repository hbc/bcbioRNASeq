#' DESeq2 Differential Expression Analysis
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2.
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @param DESeqDataSet `DESeqDataSet`.
#' @param DESeqTransform `DESeqTransform`.
#' @param DESeqResults `list`. One or more `DESeqResults`. We're using a list
#'   here to support multiple pairwise contrasts.
#' @param lfcShrink `list`. One or more `DESeqResults` returned from
#'   [DESeq2::lfcShrink()].
#'
#' @return `DESeqAnalysis`.
#'
#' @examples
#' x <- DESeqAnalysis(
#'     DESeqDataSet = dds_small,
#'     DESeqTransform = vst_small,
#'     DESeqResults = list(
#'         contrast1 = res_small
#'     ),
#'     lfcShrink = list(
#'         contrast1 = DESeq2::lfcShrink(dds = dds_small, coef = 2L)
#'     )
#' )
#' class(x)
#' slotNames(x)
DESeqAnalysis <- function(
    DESeqDataSet,
    DESeqTransform,
    DESeqResults,
    lfcShrink
) {
    new(
        Class = "DESeqAnalysis",
        DESeqDataSet = DESeqDataSet,
        DESeqTransform = DESeqTransform,
        DESeqResults = DESeqResults,
        lfcShrink = lfcShrink
    )
}
