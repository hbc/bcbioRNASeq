#' @name plotRRNAMappingRate
#' @inherit bioverbs::plotRRNAMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotRRNAMappingRate(bcb_small)
NULL



#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @importFrom bioverbs plotRRNAMappingRate
#' @usage plotRRNAMappingRate(object, ...)
#' @export
NULL



plotRRNAMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.1,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "rRNA mapping rate"
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

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("rrnaRate") * 100L,
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
                y = "rRNA mapping rate (%)",
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



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    f = "plotRRNAMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotRRNAMappingRate.bcbioRNASeq
)
