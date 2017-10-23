#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotExonicMappingRate(bcb)
#'
#' \dontrun{
#' plotExonicMappingRate(bcb, interestingGroups = "group")
#'
#' # data.frame
#' plotExonicMappingRate(metrics(bcb))
#' }
NULL



# Constructors ====
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot labs ylim
#' @importFrom viridis scale_fill_viridis
.plotExonicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 60,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~exonicRate * 100,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)") +
        ylim(0, 100)
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(fill)) {
        p <- p + scale_fill_viridis(discrete = TRUE)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotExonicMappingRate
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 60,
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotExonicMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("data.frame"),
    .plotExonicMappingRate)
