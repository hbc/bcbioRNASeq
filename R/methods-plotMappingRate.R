#' Plot Mapping Rate
#'
#' @name plotMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotMappingRate(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotMappingRate
#' @export
setMethod(
    "plotMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 90L,
        warnLimit = 70L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "mapping rate"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertIsAnImplicitInteger(passLimit)
        assert_all_are_non_negative(passLimit)
        assertIsAnImplicitInteger(warnLimit)
        assert_all_are_non_negative(warnLimit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        metrics <- metrics(object) %>%
            uniteInterestingGroups(interestingGroups)

        p <- ggplot(
            data = metrics,
            mapping = aes_(
                x = ~sampleName,
                y = ~mappedReads / totalReads * 100L,
                fill = ~interestingGroups
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            ylim(0L, 100L) +
            labs(
                title = title,
                x = "sample",
                y = "mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(passLimit)) {
            p <- p + .qcPassLine(passLimit)
        }
        if (is_positive(warnLimit)) {
            p <- p + .qcWarnLine(warnLimit)
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
