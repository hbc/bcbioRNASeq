#' Plot Intronic Mapping Rate
#'
#' @rdname plotIntronicMappingRate
#' @name plotIntronicMappingRate
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
#' plotIntronicMappingRate(bcb)
#' plotIntronicMappingRate(bcb, interestingGroups = "group")
#'
#' \dontrun{
#' # data.frame
#' metrics(bcb) %>%
#'     plotIntronicMappingRate()
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotIntronicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 20,
    flip = TRUE) {
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~intronicRate * 100,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "intronic mapping rate",
             x = "sample",
             y = "intronic mapping rate (%)") +
        ylim(0, 100) +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(warnLimit)) {
        p <- p +
            qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        warnLimit = 20,
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotIntronicMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            flip = flip)
    })



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("data.frame"),
    .plotIntronicMappingRate)
