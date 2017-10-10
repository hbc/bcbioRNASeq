#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotExonicMappingRate(bcb)
#' plotExonicMappingRate(bcb, interestingGroup = "group")
#'
#' \dontrun{
#' # data.frame
#' metrics(bcb) %>%
#'     plotExonicMappingRate()
#' }
NULL



# Constructors ====
.plotExonicMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 60,
    flip = TRUE) {
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~exonicRate * 100,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)") +
        ylim(0, 100) +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p +
            qcPassLine(passLimit)
    }
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "bcbioRNASeqANY", function(
    object,
    interestingGroup,
    passLimit = 60,
    flip = TRUE) {
    if (is.null(metrics(object))) {
        return(NULL)
    }
    if (missing(interestingGroup)) {
        interestingGroup <- interestingGroups(object)[[1]]
    }
    .plotExonicMappingRate(
        metrics(object),
        interestingGroup = interestingGroup,
        passLimit = passLimit,
        flip = flip)
})



#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "data.frame", .plotExonicMappingRate)
