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
#' plotGeneSaturation(bcb_small, label = FALSE)
#' plotGeneSaturation(bcb_small, label = TRUE)
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
        label = FALSE,
        trendline = FALSE,
        color = NULL,
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
        assert_is_a_bool(label)
        assert_is_a_bool(trendline)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsAStringOrNULL(title)

        counts <- counts(object, normalized = normalized)
        data <- metrics(object) %>%
            mutate(
                mappedReadsPerMillion = !!sym("mappedReads") / 1e6L,
                geneCount = colSums(!!counts >= !!minCounts)
            )

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "mappedReadsPerMillion",
                y = "geneCount",
                color = "interestingGroups"
            )
        ) +
            geom_point(size = 3L) +
            labs(
                title = title,
                x = "mapped reads per million",
                y = "gene count",
                color = paste(interestingGroups, collapse = ":\n")
            )

        if (isTRUE(label)) {
            p <- p + bcbio_geom_label_repel(
                mapping = aes_string(label = "sampleName")
            )
        }

        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        p
    }
)
