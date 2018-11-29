#' Plot Dispersion Estimates
#'
#' @name plotDispEsts
#' @author Michael Steinbaugh
#' @inherit DESeq2::plotDispEsts
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @details
#' This plot shows the dispersion by mean of normalized counts. We expect the
#' dispersion to decrease as the mean of normalized counts increases.
#'
#' Here we're generating a `DESeqDataSet` object on the fly, which already has
#' method support for plotting dispersion, provided by the DESeq2 package.
#'
#' @param object Object.
#'
#' @seealso `DESeq2::plotDispEsts()`.
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



#' @importFrom BiocGenerics plotDispEsts
#' @aliases NULL
#' @export
BiocGenerics::plotDispEsts



plotDispEsts.bcbioRNASeq <-  # nolint
    function() {
        validObject(object)
        # Warn and early return if any samples are duplicated.
        if (!areSamplesUnique(object)) {
            warning("Duplicate samples detected. Skipping plot.")
            return(invisible())
        }
        dds <- as(object, "DESeqDataSet")
        # Expecting warning about empty design formula.
        dds <- suppressWarnings(DESeq(dds))
        do.call(
            what = plotDispEsts,
            args = matchArgsToDoCall(args = list(object = dds))
        )
    }

formals(plotDispEsts.bcbioRNASeq) <-
    methodFormals(
        f = "plotDispEsts",
        signature = "DESeqDataSet",
        package = "DESeq2"
    )



#' @rdname plotDispEsts
#' @export
setMethod(
    f = "plotDispEsts",
    signature = signature("bcbioRNASeq"),
    definition = plotDispEsts.bcbioRNASeq
)
