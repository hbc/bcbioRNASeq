#' Plot Gene Detection Saturation
#'
#' @name plotGeneSaturation
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param trendline Include a trendline for each group.
#'
#' @return `ggplot`.
#'
#' @examples
#' # Note that minimal example distorts the y-axis
#' plotGeneSaturation(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = c("tpm", "tmm"),
        interestingGroups,
        minCounts = 1L,
        trendline = FALSE,
        color = scale_color_hue(),
        title = "gene saturation"
    ) {
        validObject(object)
        normalized <- match.arg(normalized)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assertIsAnImplicitInteger(minCounts)
        assert_all_are_in_range(minCounts, lower = 1L, upper = Inf)
        assert_is_a_bool(trendline)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        counts <- counts(object, normalized = normalized)

        p <- ggplot(
            data = metrics(object),
            mapping = aes_(
                x = ~mappedReads / 1e6L,
                y = colSums(counts >= minCounts),
                color = ~interestingGroups
            )
        ) +
            geom_point(size = 3L) +
            labs(
                title = title,
                x = "mapped reads per million",
                y = "gene count",
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        p
    }
)
