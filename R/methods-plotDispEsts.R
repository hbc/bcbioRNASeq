#' Plot Dispersion Estimates
#'
#' @rdname plotDispEsts
#' @name plotDispEsts
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#'
#' @seealso [DESeq2::plotDispEsts()].
#'
#' @return [ggplot].
#' @export
#'
#' @examples
#' data(bcb)
#' plotDispEsts(bcb)
NULL



# Methods ====
#' @rdname plotDispEsts
#' @export
setMethod("plotDispEsts", "bcbioRNADataSet", function(object) {
    bcbio(object, "DESeqDataSet") %>% plotDispEsts
})
