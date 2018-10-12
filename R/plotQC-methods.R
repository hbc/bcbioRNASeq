#' @inherit basejump::plotQC
#' @name plotQC
#' @author Michael Steinbaugh
#' @importFrom basejump plotQC
#' @export
#'
#' @examples
#' data(bcb_small)
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
    definition = .plotQC.bcbioRNASeq
)
