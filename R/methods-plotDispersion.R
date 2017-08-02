#' Plot Dispersion Estimates
#'
#' [DESeq2::plotDispEsts()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plotDispersion
#' @author Michael Steinbaugh
#' @family DESeq2 Utilities
#'
#' @param object Object.
#'
#' @seealso [DESeq2::plotDispEsts()].
#'
#' @return Dispersion plot [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plotDispersion(bcb)
setMethod("plotDispersion", "bcbioRNADataSet", function(object) {
    bcbio(object, "DESeqDataSet") %>% plotDispEsts
})
