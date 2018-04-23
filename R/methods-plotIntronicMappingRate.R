#' Plot Intronic Mapping Rate
#'
#' @name plotIntronicMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotIntronicMappingRate(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 20L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "intronic mapping rate"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        metrics <- metrics(object) %>%
            uniteInterestingGroups(interestingGroups)

        p <- ggplot(
            data = metrics,
            mapping = aes_(
                x = ~sampleName,
                y = ~intronicRate * 100L,
                fill = ~interestingGroups
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = "sample",
                y = "intronic mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            ) +
            ylim(0L, 100L)

        if (is_positive(limit)) {
            p <- p + .qcLine(limit)
        }

        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }

        if (isTRUE(flip)) {
            p <- p + coord_flip()
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = FALSE)
        }

        p
    }
)
