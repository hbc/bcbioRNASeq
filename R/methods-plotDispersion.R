#' Plot Dispersion Estimates
#'
#' [DESeq2::plotDispEsts()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plot_dispersion
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
#' plot_dispersion(bcb)
setMethod("plot_dispersion", "bcbioRNADataSet", function(object) {
    bcbio(object, "DESeqDataSet") %>% plotDispEsts
})
