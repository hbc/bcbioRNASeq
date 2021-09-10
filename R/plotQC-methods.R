#' @name plotQC
#' @author Michael Steinbaugh
#' @inherit AcidGenerics::plotQC
#' @note Updated 2021-09-03.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotQC(bcb)
NULL



## Updated 2020-09-15.
`plotQC,bcbioRNASeq` <-  # nolint
    function(object) {
        validObject(object)
        plotlist <- list(
            "totalReads" = plotTotalReads(object),
            "mappingRate" = plotMappingRate(object),
            "exonicMappingRate" = plotExonicMappingRate(object),
            "intronicMappingRate" = plotIntronicMappingRate(object),
            "rrnaMappingRate" = plotRRNAMappingRate(object),
            "x5Prime3PrimeBias" = plot5Prime3PrimeBias(object),
            "featuresDetected" = plotFeaturesDetected(object),
            "countsPerFeature" = plotCountsPerFeature(object)
        )
        if (!.isFastMode(object)) {
            plotlist[["pca"]] <- plotPCA(object)
        }
        wrap_plots(plotlist)
    }



#' @rdname plotQC
#' @export
setMethod(
    f = "plotQC",
    signature = signature("bcbioRNASeq"),
    definition = `plotQC,bcbioRNASeq`
)
