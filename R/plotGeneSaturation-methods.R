#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotGeneSaturation
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @param trendline `boolean`. Include a trendline for each group.
#'
#' @examples
#' data(bcb)
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



#' @importFrom basejump plotGeneSaturation
#' @aliases NULL
#' @export
basejump::plotGeneSaturation



plotGeneSaturation.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        minCounts = 1L,
        perMillion = TRUE,
        trendline = FALSE,
        label,
        color,
        title = "gene saturation"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsAnImplicitInteger(minCounts)
        assert_all_are_in_range(minCounts, lower = 1L, upper = Inf)
        assert_is_a_bool(perMillion)
        assert_is_a_bool(trendline)
        assert_is_a_bool(label)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsStringOrNULL(title)

        counts <- counts(object, normalized = FALSE)
        data <- metrics(object)
        assert_are_identical(colnames(counts), data[["sampleID"]])
        data[["geneCount"]] <- colSums(counts >= minCounts)

        # Convert to per million, if desired.
        xLab <- "mapped reads"
        if (isTRUE(perMillion)) {
            data <- mutate(data, mappedReads = !!sym("mappedReads") / 1e6L)
            xLab <- paste(xLab, "per million")
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("mappedReads"),
                y = !!sym("geneCount"),
                colour = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            scale_y_continuous(breaks = pretty_breaks()) +
            expand_limits(x = 0L, y = 0L) +
            labs(
                title = title,
                x = xLab,
                y = "gene count",
                colour = paste(interestingGroups, collapse = ":\n")
            )

        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(label)) {
            p <- p + basejump_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }

        p
    }

formals(plotGeneSaturation.bcbioRNASeq)[["color"]] <-
    formalsList[["color.discrete"]]
formals(plotGeneSaturation.bcbioRNASeq)[["label"]] <-
    formalsList[["label"]]



#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature("bcbioRNASeq"),
    definition = plotGeneSaturation.bcbioRNASeq
)
