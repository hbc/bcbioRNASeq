#' Plot Dispersion Estimates
#'
#' This plot shows the dispersion by mean of normalized counts. We expect the
#' dispersion to decrease as the mean of normalized counts increases.
#'
#' @note Here we're generating a `DESeqDataSet` object on the fly, which already
#'   has method support for plotting dispersion, provided by the DESeq2 package.
#'
#' @name plotDispEsts
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @importFrom BiocGenerics plotDispEsts
#' @export
#'
#' @inherit DESeq2::plotDispEsts
#' @param object Object.
#'
#' @seealso
#' - [DESeq2::plotDispEsts()].
#' - `getMethod("plotDispEsts", "DESeqDataSet")`.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotDispEsts(bcb_small)
#'
#' # Custom colors, using DESeq2 parameters.
#' plotDispEsts(
#'     object = bcb_small,
#'     genecol = "gray",
#'     fitcol = "purple",
#'     finalcol = "orange"
#' )
NULL



.plotDispEsts.bcbioRNASeq <-  # nolint
    function() {
        validObject(object)
        dds <- as(object, "DESeqDataSet")
        dds <- suppressWarnings(DESeq(dds))
        args <- setArgsToDoCall(
            args = list(object = dds),
            call = matchCall()
        )
        do.call(what = plotDispEsts, args = args)
    }
formals(.plotDispEsts.bcbioRNASeq) <- methodFormals(
    f = "plotDispEsts",
    signature = "DESeqDataSet"
)



#' @rdname plotDispEsts
#' @export
setMethod(
    f = "plotDispEsts",
    signature = signature("bcbioRNASeq"),
    definition = .plotDispEsts.bcbioRNASeq
)
