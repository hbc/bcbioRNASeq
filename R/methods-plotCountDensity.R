#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @name plotCountDensity
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#' @param style Desired plot style (`color` or `fill`).
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotCountDensity(bcb)
#' plotCountDensity(bcb, interestingGroup = "group")
#'
#' \dontrun{
#' # data.frame
#' meltLog10(bcb, normalized = "tmm") %>%
#'     plotCountDensity()
#' }
NULL



# Constructors ====
.plotCountDensity <- function(
    object,
    interestingGroup = "sampleName",
    style = "color") {
    validStyles <- c("color", "fill")
    if (!style %in% validStyles) {
        stop(paste(
            "Valid 'style' arguments:",
            toString(validStyles)))
    }
    if (style == "color") {
        color <- interestingGroup
        fill <- NULL
    } else if (style == "fill") {
        color <- NULL
        fill <- interestingGroup
    }
    p <- ggplot(
        object,
        mapping = aes_string(
            x = "counts",
            group = interestingGroup,
            color = color,
            fill = fill)
    )
    if (style == "color") {
        p <- p +
            geom_density()
    } else if (style == "fill") {
        p <- p +
            geom_density(alpha = 0.75, color = NA)
    }
    p +
        labs(title = "count density",
             x = "log10 counts per gene") +
        scale_color_viridis(discrete = TRUE) +
        scale_fill_viridis(discrete = TRUE)
}



# Methods ====
#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroup,
        normalized = "tmm",
        style = "color") {
        if (missing(interestingGroup)) {
            interestingGroup <- interestingGroups(object)[[1]]
        }
        .plotCountDensity(
            meltLog10(object, normalized = normalized),
            interestingGroup = interestingGroup,
            style = style)
    })



#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("data.frame"),
    .plotCountDensity)
