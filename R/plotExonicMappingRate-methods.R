#' @name plotExonicMappingRate
#' @inherit bioverbs::plotExonicMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @description Ideally, at least 60 percent of total reads should map to exons.
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotExonicMappingRate(bcb_small)
NULL



#' @rdname plotExonicMappingRate
#' @name plotExonicMappingRate
#' @importFrom bioverbs plotExonicMappingRate
#' @usage plotExonicMappingRate(object, ...)
#' @export
NULL



plotExonicMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups,
        limit = 0.6,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "exonic mapping rate"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assert_is_a_number(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        p <- metrics(object) %>%
            ggplot(
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("exonicRate") * 100L,
                    fill = !!sym("interestingGroups")
                )
            ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = "exonic mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            # Convert to percentage
            if (limit > 1L) {
                # nocov start
                warning("`limit`: Use ratio (0-1) instead of percentage")
                # nocov end
            } else {
                limit <- limit * 100L
            }
            if (limit < 100L) {
                p <- p + bcbio_geom_abline(yintercept = limit)
            }
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



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    f = "plotExonicMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotExonicMappingRate.bcbioRNASeq
)
