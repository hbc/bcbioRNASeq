#' Plot Intronic Mapping Rate
#'
#' @rdname plotIntronicMappingRate
#' @name plotIntronicMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotIntronicMappingRate(bcb)
#'
#' \dontrun{
#' plotIntronicMappingRate(
#'     bcb,
#'     interestingGroups = "group"
#'     fill = NULL)
#' }
#'
#' # data.frame
#' \dontrun{
#' metrics(bcb) %>% plotIntronicMappingRate()
#' }
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom viridis scale_fill_viridis
.plotIntronicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 20,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~intronicRate * 100,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(title = "intronic mapping rate",
             x = "sample",
             y = "intronic mapping rate (%)",
             fill = paste(interestingGroups, collapse = ":\n")) +
        ylim(0, 100)
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (!is.null(fill)) {
        p <- p + fill
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotIntronicMappingRate
#' @importFrom viridis scale_color_viridis
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 20,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]]
        }
        .plotIntronicMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("data.frame"),
    .plotIntronicMappingRate)
