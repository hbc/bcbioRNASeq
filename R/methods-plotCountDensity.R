#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @name plotCountDensity
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNADataSet
#' plotCountDensity(bcb)
#'
#' # data.frame
#' meltLog10(bcb) %>% plotCountDensity
NULL



# Constructors ====
.plotCountDensity <- function(
    object,
    interestingGroup = "sampleName") {
    ggplot(object,
        aes_(x = ~counts,
             group = as.name(interestingGroup),
             color = as.name(interestingGroup))) +
        geom_density() +
        labs(title = "count density",
             x = "log10 counts per gene") +
        scale_color_viridis(discrete = TRUE)
}



# Methods ====
#' @rdname plotCountDensity
#' @export
setMethod("plotCountDensity", "bcbioRNADataSet", function(
    object, normalized = "tmm") {
    .plotCountDensity(
        meltLog10(object, normalized = normalized),
        interestingGroup = .interestingGroup(object))
})



#' @rdname plotCountDensity
#' @export
setMethod("plotCountDensity", "data.frame", .plotCountDensity)
