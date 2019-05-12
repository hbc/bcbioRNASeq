#' @name plotGeneSaturation
#' @inherit bioverbs::plotGeneSaturation
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param trendline `boolean`. Include a trendline for each group.
#' @param ... Additional arguments.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotGeneSaturation(bcb_small, label = FALSE)
#' plotGeneSaturation(bcb_small, label = TRUE)
NULL



#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @importFrom bioverbs plotGeneSaturation
#' @usage plotGeneSaturation(object, ...)
#' @export
NULL



plotGeneSaturation.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups,
        minCounts = 1L,
        trendline = FALSE,
        label = getOption("bcbio.label", FALSE),
        color = getOption("bcbio.discrete.color", NULL),
        title = "gene saturation"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assertIsAnImplicitInteger(minCounts)
        assert_all_are_in_range(minCounts, lower = 1L, upper = Inf)
        assert_is_a_bool(trendline)
        assert_is_a_bool(label)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        counts <- counts(object, normalized = FALSE)
        p <- metrics(object) %>%
            mutate(geneCount = colSums(!!counts >= !!minCounts)) %>%
            ggplot(
                mapping = aes(
                    x = !!sym("mappedReads") / 1e6L,
                    y = !!sym("geneCount"),
                    color = !!sym("interestingGroups")
                )
            ) +
            geom_point(size = 3L) +
            scale_y_continuous(breaks = pretty_breaks()) +
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

        if (isTRUE(label)) {
            p <- p + bcbio_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }

        p
    }



#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature("bcbioRNASeq"),
    definition = plotGeneSaturation.bcbioRNASeq
)
