#' @name plotQc
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::plotQc
#' @note Updated 2022-03-07.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotQc(bcb)
NULL



## Updated 2020-09-15.
`plotQc,bcbioRNASeq` <- # nolint
    function(object) {
        validObject(object)
        plotlist <- list(
            "totalReads" = plotTotalReads(object),
            "mappingRate" = plotMappingRate(object),
            "exonicMappingRate" = plotExonicMappingRate(object),
            "intronicMappingRate" = plotIntronicMappingRate(object),
            "rrnaMappingRate" = plotRrnaMappingRate(object),
            "x5Prime3PrimeBias" = plot5Prime3PrimeBias(object),
            "featuresDetected" = plotFeaturesDetected(object),
            "countsPerFeature" = plotCountsPerFeature(object)
        )
        if (!.isFastMode(object)) {
            plotlist[["pca"]] <- plotPca(object)
        }
        wrap_plots(plotlist)
    }



#' @rdname plotQc
#' @export
setMethod(
    f = "plotQc",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotQc,bcbioRNASeq`
)
