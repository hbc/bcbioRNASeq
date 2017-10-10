#' Plot Mapping Rate
#'
#' @rdname plotMappingRate
#' @name plotMappingRate
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
#' plotMappingRate(bcb)
#'
#' \dontrun{
#' plotMappingRate(bcb, interestingGroups = "group")
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plotMappingRate()
#' }
NULL



# Constructors ====
.plotMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    passLimit = 90,
    warnLimit = 70,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~mappedReads / totalReads * 100,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        ylim(0, 100) +
        labs(title = "mapping rate",
             x = "sample",
             y = "mapping rate (%)") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("bcbioRNASeqANY"),
    function(
        object,
        passLimit = 90,
        warnLimit = 70,
        flip = TRUE) {
        .plotMappingRate(
            metrics(object),
            interestingGroups = interestingGroups(object)[[1]],
            passLimit = passLimit,
            warnLimit = warnLimit,
            flip = flip)
    })



#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("data.frame"),
    .plotMappingRate)
