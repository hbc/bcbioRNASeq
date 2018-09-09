#' DESeq2 Differential Expression Analysis
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2.
#'
#' @section DESeqDataSet:
#'
#' We recommend generating the `DESeqDataSet` by coercion from `bcbioRNASeq`
#' object using `as(dds, "bcbioRNASeq")`. Don't use the [DESeq2::DESeqDataSet()]
#' or [DESeq2::DESeqDataSetFromMatrix()] constructors to generate the
#' `DESeqDataSet` object.
#'
#' @section DESeqResults:
#'
#' Don't modify any of the `DESeqResults` objects manually. This includes
#' rearranging the rows or dropping genes without adjusted P values. We'll take
#' care of this automatically in supported methods.
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
#' library(DESeq2)
#' dds <- dds_small
#' class(dds)
#' vst <- varianceStabilizingTransformation(dds)
#' class(vst)
#' resultsNames(dds)
#' res <- results(dds, name = resultsNames(dds)[[2L]])
#' class(res)
#' x <- DESeqAnalysis(
#'     DESeqDataSet = dds_small,
#'     DESeqTransform = vst_small,
#'     DESeqResults = list(
#'         res_small
#'     ),
#'     lfcShrink = list(
#'         lfcShrink(dds = dds_small, coef = 2L)
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
