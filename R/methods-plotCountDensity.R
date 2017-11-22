#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @name plotCountDensity
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotGene
#'
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param style Desired plot style (`line` or `solid`).
#'
#' @examples
#' # bcbioRNASeq
#' bcb <- examples[["bcb"]]
#' plotCountDensity(bcb)
#' \dontrun{
#' plotCountDensity(
#'     bcb,
#'     interestingGroups = "group",
#'     fill = NULL)
#' }
#'
#' # data.frame
#' \dontrun{
#' meltLog10(bcb, normalized = "tmm") %>%
#'     plotCountDensity()
#' }
NULL



# Constructors ====
#' @importFrom basejump uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_density ggplot labs
#' @importFrom viridis scale_color_viridis scale_fill_viridis
.plotCountDensity <- function(
    object,
    interestingGroups = "sampleName",
    style = "solid",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE)) {
    validStyles <- c("line", "solid")
    if (!style %in% validStyles) {
        stop(paste(
            "Valid 'style' arguments:",
            toString(validStyles)
        ), call. = FALSE)
    }
    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_string(
            x = "counts",
            group = "interestingGroups",
            color = "interestingGroups",
            fill = "interestingGroups")
    ) +
        labs(title = "count density",
             x = "log10 counts per gene",
             fill = paste(interestingGroups, collapse = ":\n"))
    if (style == "line") {
        p <- p +
            geom_density(fill = NA) +
            color
    } else if (style == "solid") {
        p <- p +
            geom_density(alpha = 0.75, color = NA) +
            fill
    }
    p
}



# Methods ====
#' @rdname plotCountDensity
#' @importFrom S4Vectors metadata
#' @importFrom viridis scale_color_viridis scale_fill_viridis
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "tmm",
        style = "solid",
        color = scale_color_viridis(discrete = TRUE),
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
             interestingGroups <- basejump::interestingGroups(object)
        }
        .plotCountDensity(
            meltLog10(object, normalized = normalized),
            interestingGroups = interestingGroups,
            style = style,
            color = color,
            fill = fill)
    })



#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("data.frame"),
    .plotCountDensity)
