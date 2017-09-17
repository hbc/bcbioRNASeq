#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
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
#' # bcbioRNADataSet
#' plotRRNAMappingRate(bcb)
#' plotRRNAMappingRate(bcb, interestingGroup = "group")
#'
#' # data.frame
#' metrics(bcb) %>% plotRRNAMappingRate
NULL



# Constructors ====
.plotRRNAMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    warnLimit = 10L,
    flip = TRUE) {
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~rRnaRate * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "rrna mapping rate",
             x = "sample",
             y = "rRNA mapping rate (%)") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotRRNAMappingRate
#' @export
setMethod("plotRRNAMappingRate", "bcbioRNADataSet", function(
    object,
    interestingGroup,
    warnLimit = 10L,
    flip = TRUE) {
    if (is.null(metrics(object))) {
        return(NULL)
    }
    if (missing(interestingGroup)) {
        interestingGroup <- .interestingGroup(object)
    }
    .plotRRNAMappingRate(
        metrics(object),
        interestingGroup = interestingGroup,
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plotRRNAMappingRate
#' @export
setMethod("plotRRNAMappingRate", "data.frame", .plotRRNAMappingRate)
