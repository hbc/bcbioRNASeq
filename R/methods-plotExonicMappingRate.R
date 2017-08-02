#' Plot Exonic Mapping Rate
#'
#' @rdname plotExonicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @family Quality Control Plots
#' @inherit qcPlots
#'
#' @examples
#' data(bcb)
#' plotExonicMappingRate(bcb)



#' @rdname plotExonicMappingRate
.plotExonicMappingRate <- function(
    metrics,
    interestingGroup = "sampleName",
    passLimit = 60L,
    flip = TRUE) {
    if (is.null(metrics)) return(NULL)
    p <- ggplot(metrics,
                aes_(x = ~sampleName,
                     y = ~exonicRate * 100L,
                     fill = as.name(interestingGroup))) +
        geom_bar(stat = "identity") +
        labs(title = "exonic mapping rate",
             x = "sample",
             y = "exonic mapping rate (%)") +
        ylim(0L, 100L)
    if (!is.null(passLimit)) {
        p <- p +
            geom_hline(alpha = qcLineAlpha,
                       color = qcPassColor,
                       size = qcLineSize,
                       yintercept = passLimit)
    }
    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }
    p
}



#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "bcbioRNADataSet", function(object, ...) {
    .plotExonicMappingRate(
        metrics = metrics(object),
        interestingGroup = .interestingGroup(object),
        ...)
})



#' @rdname plotExonicMappingRate
#' @export
setMethod("plotExonicMappingRate", "data.frame", .plotExonicMappingRate)
