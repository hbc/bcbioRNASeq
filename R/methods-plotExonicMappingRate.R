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
#'
#' \dontrun{
#' plotExonicMappingRate(bcb, interestingGroups = "group")
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plotExonicMappingRate()
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotExonicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 60,
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
setMethod(
    "plotExonicMappingRate",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        passLimit = 60,
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
            flip = flip)
    })



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("data.frame"),
    .plotExonicMappingRate)
