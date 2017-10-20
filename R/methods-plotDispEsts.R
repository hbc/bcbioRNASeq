#' Plot Dispersion Estimates
#'
#' @rdname plotDispEsts
#' @name plotDispEsts
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotDispEsts
#'
#' @inheritParams AllGenerics
#' @param ... Passthrough arguments to [DESeq2::plotDispEsts()].
#'
#' @seealso [DESeq2::plotDispEsts()].
#'
#' @return [ggplot].
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotDispEsts(bcb)
#'
#' \dontrun{
#' plotDispEsts(
#'     bcb,
#'     genecol = "gray",
#'     fitcol = "purple",
#'     finalcol = "orange")
#' }
NULL



# Methods ====
#' @rdname plotDispEsts
#' @export
setMethod(
    "plotDispEsts",
    signature("bcbioRNASeqANY"),
    function(object, ...) {
        bcbio(object, "DESeqDataSet") %>%
            plotDispEsts(...)
    })
