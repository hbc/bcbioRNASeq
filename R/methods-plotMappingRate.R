#' Plot Mapping Rate
#'
#' @rdname plotMappingRate
#' @name plotMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' plotMappingRate(bcb)
#'
#' \dontrun{
#' plotMappingRate(
#'     bcb,
#'     interestingGroups = "group",
#'     fill = NULL)
#' }
#'
#' # data.frame
#' \dontrun{
#' metrics(bcb) %>% plotMappingRate()
#' }
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom viridis scale_fill_viridis
.plotMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 90,
    warnLimit = 70,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE) {
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / totalReads * 100,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        ylim(0, 100) +
        labs(title = "mapping rate",
             x = "sample",
             y = "mapping rate (%)",
             fill = paste(interestingGroups, collapse = ":\n"))
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
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
#' @rdname plotMappingRate
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 90,
        warnLimit = 70,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        .plotMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            passLimit = passLimit,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip)
    })



#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("data.frame"),
    .plotMappingRate)
