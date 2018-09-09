#' DESeq2 Differential Expression Analysis
#'
#' Class containing all elements generated during differential expression
#' analysis with DESeq2.
#'
#' @family S4 Object
#' @author Michael Steinbaugh
#' @export
#'
#' @param data `DESeqDataSet`.
#' @param transform `DESeqTransform`.
#' @param results `list`. One or more `DESeqResults`. We're using a list here
#'   to support multiple pairwise contrasts.
#'
#' @examples
#' x <- DESeqAnalysis(
#'     data = dds_small,
#'     transform = vst_small,
#'     results = list(res_small)
#' )
#' class(x)
#' slotNames(x)
DESeqAnalysis <- function(
    data,
    transform,
    results
) {
    new(
        Class = "DESeqAnalysis",
        data = data,
        transform = transform,
        results = results
    )
}
