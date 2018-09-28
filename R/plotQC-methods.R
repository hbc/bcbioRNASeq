#' Plot Quality Control
#'
#' @name plotQC
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @importFrom basejump plotQC
#' @inherit basejump::plotQC
#' @export
#'
#' @examples
#' plotQC(bcb_small)
NULL



.plotQC.bcbioRNASeq <-  # nolint
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
                plotGeneSaturation(object),
                plotCountsPerGene(object),
                plotPCA(object)
            )
        )
    }
