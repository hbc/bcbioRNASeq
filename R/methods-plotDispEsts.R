#' Plot Dispersion Estimates
#'
#' Method support for plotting the dispersion of counts stored in a
#' [bcbioRNASeq] object. Here we're using the internally stored [DESeqDataSet],
#' which already has method support for plotting dispersion, provided by
#' the DESeq2 package.
#'
#' @rdname plotDispEsts
#' @name plotDispEsts
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotDispEsts
#'
#' @param object [bcbioRNASeq].
#' @param ... Passthrough arguments to [DESeq2::plotDispEsts()].
#'
#' @seealso
#' - [DESeq2::plotDispEsts()].
#' - `getMethod("plotDispEsts", "DESeqDataSet")`.
#'
#' @return [ggplot].
#'
#' @examples
#' bcb <- examples[["bcb"]]
#' plotDispEsts(bcb)
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
    signature("bcbioRNASeq"),
    function(object, ...) {
        bcbio(object, "DESeqDataSet") %>%
            plotDispEsts(...)
    })
