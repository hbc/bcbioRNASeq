#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @name plotCountDensity
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotCountDensity(bcb, normalized = "tmm")
#'
#' # data.frame
#' meltLog10(bcb, normalized = "tmm") %>%
#'     plotCountDensity
NULL



# Constructors ====
.plotCountDensity <- function(
    object,
    interestingGroup = "sampleName") {
    ggplot(object,
        aes_(x = ~counts,
             group = as.name(interestingGroup),
             fill = as.name(interestingGroup))) +
        geom_density(alpha = 0.75, color = NA) +
        labs(title = "count density",
             x = "log10 counts per gene") +
        scale_fill_viridis(discrete = TRUE)
}



# Methods ====
#' @rdname plotCountDensity
#' @export
setMethod("plotCountDensity", "bcbioRNADataSet", function(
    object,
    normalized = "tmm") {
    .plotCountDensity(
        meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object))
})



#' @rdname plotCountDensity
#' @export
setMethod("plotCountDensity", "data.frame", .plotCountDensity)
