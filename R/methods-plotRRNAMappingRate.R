#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotRRNAMappingRate(bcb)
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
    if (is.null(object)) return(NULL)
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
    warnLimit = 10L,
    flip = TRUE) {
    .plotRRNAMappingRate(
        metrics(object),
        interestingGroup = .interestingGroup(object),
        warnLimit = warnLimit,
        flip = flip)
})



#' @rdname plotRRNAMappingRate
#' @export
setMethod("plotRRNAMappingRate", "data.frame", .plotRRNAMappingRate)
