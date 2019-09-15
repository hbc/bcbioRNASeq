#' @name plotCountsPerFeature
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit acidplots::plotCountsPerFeature
#' @note Updated 2019-09-15.
#'
#' @inheritParams plotCounts
#' @inheritParams acidroxygen::params
#' @param ... Passthrough to `SummarizedExperiment` method defined in acidplots.
#'   See [acidplots::plotCountsPerFeature()] for details.
#'
#' @section Trimmed Mean of M-Values:
#'
#' We recommend visualizing counts normalized with the **T**rimmed **M**ean of
#' **M**-Values (TMM) method here. TMM normalization equates the overall
#' expression levels of genes between samples under the assumption that the
#' majority of them are not differentially expressed. Therefore, by normalizing
#' for total RNA expression by sample, we expect the spread of the
#' TMM-normalized counts per gene to be similar for every sample.
#'
#' @references TMM: Robinson, et al., 2010.
#'
#' @examples
#' data(bcb)
#' plotCountsPerFeature(bcb)
NULL



#' @rdname plotCountsPerFeature
#' @name plotCountsPerFeature
#' @importFrom bioverbs plotCountsPerFeature
#' @importMethodsFrom acidplots plotCountsPerFeature
#' @usage plotCountsPerFeature(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotCountsPerFeature,bcbioRNASeq` <-  # nolint
    function(object, normalized, ...) {
        normalized <- match.arg(normalized)
        args <- .dynamicTrans(object = object, normalized = normalized, ...)
        do.call(what = plotCountsPerFeature, args = args)
    }

f <- formals(`plotCountsPerFeature,bcbioRNASeq`)
f[["normalized"]] <- unique(c("tmm", normalizedCounts))
formals(`plotCountsPerFeature,bcbioRNASeq`) <- f



#' @rdname plotCountsPerFeature
#' @export
setMethod(
    f = "plotCountsPerFeature",
    signature = signature("bcbioRNASeq"),
    definition = `plotCountsPerFeature,bcbioRNASeq`
)



## Note that this function is defined in F1000 v2 workflow paper.

#' @rdname plotCountsPerFeature
#' @export
plotCountDensity <- function(object, ...) {
    plotCountsPerFeature(object = object, geom = "density", ...)
}
