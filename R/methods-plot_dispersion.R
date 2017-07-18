#' Plot dispersion estimates
#'
#' [DESeq2::plotDispEsts()] wrapper supporting a [bcbioRNADataSet].
#'
#' @rdname plot_dispersion
#' @family DESeq2 utilities
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @examples
#' data(bcb)
#' plot_dispersion(bcb)
#' @return Dispersion plot [ggplot].
#' @export
setMethod("plot_dispersion", "bcbioRNADataSet", function(object) {
    bcbio(object, "DESeqDataSet") %>% plotDispEsts
})
