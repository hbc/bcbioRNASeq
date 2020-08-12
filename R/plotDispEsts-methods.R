#' Plot dispersion estimates
#'
#' @details
#' This plot shows the dispersion by mean of normalized counts. We expect the
#' dispersion to decrease as the mean of normalized counts increases.
#'
#' Here we're generating a `DESeqDataSet` object on the fly, which already has
#' method support for plotting dispersion, provided by the DESeq2 package.
#'
#' @name plotDispEsts
#' @author Michael Steinbaugh
#' @inherit DESeq2::plotDispEsts description
#' @note Updated 2019-09-15.
#'
#' @param object Object.
#' @param ... Passthrough to `DESeqDataSet` method defined in DESeq2.
#'   See [DESeq2::plotDispEsts()] for details.
#'
#' @seealso [DESeq2::plotDispEsts()].
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb)
#' plotDispEsts(bcb)
#'
#' ## Custom colors, using DESeq2 parameters.
#' plotDispEsts(
#'     object = bcb,
#'     genecol = "gray",
#'     fitcol = "purple",
#'     finalcol = "orange"
#' )
NULL



#' @rdname plotDispEsts
#' @name plotDispEsts
#' @importFrom BiocGenerics plotDispEsts
#' @usage plotDispEsts(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotDispEsts,bcbioRNASeq` <-  # nolint
    function(object, ...) {
        validObject(object)
        ## Warn and early return if any samples are duplicated.
        if (!hasUniqueCols(object)) {
            warning("Duplicate samples detected. Skipping plot.")
            return()
        }
        dds <- as(object, "DESeqDataSet")
        ## Expecting warning about empty design formula.
        suppressWarnings(dds <- DESeq(dds))
        plotDispEsts(object = dds, ...)
    }



#' @rdname plotDispEsts
#' @export
setMethod(
    f = "plotDispEsts",
    signature = signature("bcbioRNASeq"),
    definition = `plotDispEsts,bcbioRNASeq`
)
