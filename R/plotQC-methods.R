#' @name plotQC
#' @author Michael Steinbaugh
#' @inherit bioverbs::plotQC
#' @examples
#' data(bcb)
#' plotQC(bcb)
NULL



#' @importFrom bioverbs plotQC
#' @aliases NULL
#' @export
bioverbs::plotQC



plotQC.bcbioRNASeq <-  # nolint
    function(object) {
        validObject(object)
        plot_grid(
            plotlist = list(
                plotTotalReads(object),
                plotMappingRate(object),
                plotExonicMappingRate(object),
                plotIntronicMappingRate(object),
                plotRRNAMappingRate(object),
                plot5Prime3PrimeBias(object),
                plotGenesDetected(object),
                plotCountsPerGene(object),
                plotPCA(object)
            ),
            nrow = 3L,
            ncol = 3L
        )
    }



#' @rdname plotQC
#' @export
setMethod(
    f = "plotQC",
    signature = signature("bcbioRNASeq"),
    definition = plotQC.bcbioRNASeq
)
