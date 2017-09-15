#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotExonicMappingRate(bcb)
#'
#' # data.frame
#' metrics(bcb) %>% plotExonicMappingRate
NULL



# Constructors ====
.plotExonicMappingRate <- function(
    object,
    interestingGroup = "sampleName",
    passLimit = 60L,
    flip = TRUE) {
    if (is.null(object)) return(NULL)
    p <- ggplot(object,
                aes_(x = ~sampleName,
                     y = ~exonicRate * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)") +
        ylim(0L, 100L) +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(passLimit)) {
        p <- p + qcPassLine(passLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "bcbioRNADataSet", function(
    object,
    passLimit = 60L,
    flip = TRUE) {
    .plotExonicMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        passLimit = passLimit,
        flip = flip)
})



#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "data.frame", .plotExonicMappingRate)
