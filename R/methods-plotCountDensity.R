#' Plot Count Density
#'
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
#' plotCountDensity(bcb_small)
#' plotCountDensity(
#'     object = bcb_small,
#'     style = "line",
#'     interestingGroups = "sampleName"
#' )
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_density ggplot guides labs
.plotCountDensity <- function(
    object,
    interestingGroups = "sampleName",
    style = "solid",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE),
    title = TRUE,
    xlab = "counts"
) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assert_is_a_string(style)
    assert_is_subset(style, c("line", "solid"))
    assertIsColorScaleDiscreteOrNULL(color)
    assertIsFillScaleDiscreteOrNULL(fill)

    # Title
    if (isTRUE(title)) {
        title <- "count density"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "counts",
            group = "interestingGroups",
            color = "interestingGroups",
            fill = "interestingGroups"
        )
    ) +
        labs(
            title = title,
            x = xlab,
            fill = paste(interestingGroups, collapse = ":\n")
        )

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
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "rlog",
        style = "solid",
        color = scale_color_viridis(discrete = TRUE),
        fill = scale_fill_viridis(discrete = TRUE),
        title = TRUE
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assert_is_a_string(normalized)
        if (normalized %in% c("rlog", "vst")) {
            fxn <- .meltCounts
        } else {
            fxn <- .meltLog2Counts
        }
        melt <- fxn(
            counts(object, normalized = normalized),
            colData = colData(object)
        )
        .plotCountDensity(
            object = melt,
            interestingGroups = interestingGroups,
            style = style,
            color = color,
            fill = fill,
            title = title,
            xlab = paste("log2", normalized, "counts")
        )
    }
)
