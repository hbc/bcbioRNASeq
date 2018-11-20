#' @name plotQC
#' @author Michael Steinbaugh
#' @inherit basejump::plotQC
#' @examples
#' data(bcb)
#' plotQC(bcb)
NULL



#' @importFrom basejump plotQC
#' @aliases NULL
#' @export
basejump::plotQC



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
