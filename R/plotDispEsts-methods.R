#' Plot dispersion estimates
#'
#' @name plotDispEsts
#' @author Michael Steinbaugh
#' @note Updated 2022-03-07.
#'
#' @details
#' This plot shows the dispersion by mean of normalized counts. We expect the
#' dispersion to decrease as the mean of normalized counts increases.
#'
#' Here we're generating a `DESeqDataSet` object on the fly, which already has
#' method support for plotting dispersion, provided by the DESeq2 package.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Passthrough to `DESeqDataSet` method defined in DESeq2.
#' See `DESeq2::plotDispEsts()` for details.
#'
#' @return `ggplot`.
#'
#' @seealso
#' - `DESeq2::plotDispEsts()`.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
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



## Updated 2019-09-15.
`plotDispEsts,bcbioRNASeq` <- # nolint
    function(object, ...) {
        validObject(object)
        ## Warn and early return if any samples are duplicated.
        if (!hasUniqueCols(object)) {
            alertWarning("Duplicate samples detected. Skipping plot.")
            return(invisible(NULL))
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
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotDispEsts,bcbioRNASeq`
)
