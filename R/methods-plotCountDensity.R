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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotCountDensity(bcb, style = "solid")
#' plotCountDensity(
#'     bcb,
#'     style = "line",
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- meltLog10(bcb, normalized = "tmm")
#' plotCountDensity(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_density ggplot guides labs
#' @importFrom viridis scale_color_viridis scale_fill_viridis
.plotCountDensity <- function(
    object,
    interestingGroups = "sampleName",
    style = "solid",
    color = viridis::scale_color_viridis(discrete = TRUE),
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    title = TRUE) {
    validStyles <- c("line", "solid")
    if (!style %in% validStyles) {
        abort(paste(
            "Valid `style` arguments:",
            toString(validStyles)
        ))
    }

    if (isTRUE(title)) {
        title <- "count density"
    } else if (!is.character(title)) {
        title <- NULL
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
        labs(title = title,
             x = "log10 counts per gene",
             fill = paste(interestingGroups, collapse = ":\n"))

    if (style == "line") {
        p <- p + geom_density(fill = NA)
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
    } else if (style == "solid") {
        p <- p + geom_density(alpha = 0.75, color = NA)
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plotCountDensity
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
        color = viridis::scale_color_viridis(discrete = TRUE),
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        title = TRUE) {
        if (missing(interestingGroups)) {
             interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotCountDensity(
            meltLog10(object, normalized = normalized),
            interestingGroups = interestingGroups,
            style = style,
            color = color,
            fill = fill,
            title = title)
    })



#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("data.frame"),
    .plotCountDensity)
