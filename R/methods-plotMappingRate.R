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
#' # data.frame
#' metrics(bcb) %>% plotMappingRate
NULL



# Constructors ====
.plotMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 90L,
    warnLimit = 70L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~mappedReads / totalReads * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        ylim(0L, 100L) +
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
setMethod("plotMappingRate", "bcbioRNASeqANY", function(
    object,
    passLimit = 90L,
    warnLimit = 70L,
    flip = TRUE) {
    .plotMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        passLimit = passLimit,
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plotMappingRate
#' @export
setMethod("plotMappingRate", "data.frame", .plotMappingRate)
