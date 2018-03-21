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
#' # Minimal exampe distorts the y-axis
#' plotGeneSaturation(bcb_small, interestingGroups = "sampleName")
NULL



# Constructors =================================================================
#' @importFrom bcbioBase interestingGroups uniteInterestingGroups
#' @importFrom ggplot2 aes_ geom_point geom_smooth ggplot labs
.plotGeneSaturation <- function(
    object,
    normalized = c("rlog", "vst", "tmm", "tpm"),
    interestingGroups,
    minCounts = 0L,
    trendline = TRUE,
    color = scale_color_viridis(discrete = TRUE),
    title = "gene saturation"
) {
    validObject(object)
    normalized <- match.arg(normalized)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assertIsAnImplicitInteger(minCounts)
    assert_all_are_non_negative(minCounts)
    assert_is_a_bool(trendline)
    assertIsColorScaleDiscreteOrNULL(color)
    assertIsAStringOrNULL(title)

    metrics <- metrics(object) %>%
        uniteInterestingGroups(interestingGroups)
    counts <- counts(object, normalized = normalized)

    p <- ggplot(
        data = metrics,
        mapping = aes_(
            x = ~mappedReads / 1e6L,
            y = colSums(counts > minCounts),
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



# Methods ======================================================================
#' @rdname plotGeneSaturation
#' @export
setMethod(
    "plotGeneSaturation",
    signature("bcbioRNASeq"),
    .plotGeneSaturation
)
