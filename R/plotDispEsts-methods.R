#' @importFrom BiocGenerics plotDispEsts
#' @aliases NULL
#' @export
BiocGenerics::plotDispEsts



#' Plot Dispersion Estimates
#' @name plotDispEsts
#' @inherit DESeq2::plotDispEsts
#' @author Michael Steinbaugh
#'
#' @details
#' This plot shows the dispersion by mean of normalized counts. We expect the
#' dispersion to decrease as the mean of normalized counts increases.
#'
#' Here we're generating a `DESeqDataSet` object on the fly, which already has
#' method support for plotting dispersion, provided by the DESeq2 package.
#'
#' @inheritParams general
#' @param object Object.
#'
#' @seealso [DESeq2::plotDispEsts].
#'
#' @return `ggplot`.
#'
#' @examples
#' data(bcb_small)
#' plotDispEsts(bcb_small)
#'
#' ## Custom colors, using DESeq2 parameters.
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
        # Expecting warning about empty design formula.
        dds <- suppressWarnings(DESeq(dds))
        do.call(
            what = plotDispEsts,
            args = matchArgsToDoCall(args = list(object = dds))
        )
    }
formals(.plotDispEsts.bcbioRNASeq) <-
    methodFormals(f = "plotDispEsts", signature = "DESeqDataSet")



#' @describeIn plotDispEsts Coerces to `DESeqDataSet` and inherits method from
#'   DESeq2.
#' @export
setMethod(
    f = "plotDispEsts",
    signature = signature("bcbioRNASeq"),
    definition = .plotDispEsts.bcbioRNASeq
)



# DESeqAnalysis ================================================================
.plotDispEsts.DESeqAnalysis <-  # nolint
    function(object, ...) {
        plotDispEsts(as(object, "DESeqDataSet"), ...)
    }
formals(.plotDispEsts.DESeqAnalysis) <-
    methodFormals(f = "plotDispEsts", signature = "DESeqDataSet")



#' @describeIn plotDispEsts Extracts internal `DESeqDataSet` and inherits method
#'   from DESeq2.
#' @export
setMethod(
    f = "plotDispEsts",
    signature = signature("DESeqAnalysis"),
    definition = .plotDispEsts.DESeqAnalysis
)
