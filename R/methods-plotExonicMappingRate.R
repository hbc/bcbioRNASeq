#' Plot Exonic Mapping Rate
#'
#' @name plotExonicMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotExonicMappingRate(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotExonicMappingRate
#' @export
setMethod(
    "plotExonicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 60L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "exonic mapping rate"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        p <- ggplot(
            data = metrics(object),
            mapping = aes_(
                x = ~sampleName,
                y = ~exonicRate * 100L,
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
                y = "exonic mapping rate (%)",
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
