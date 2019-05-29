#' @name plotTotalReads
#' @inherit bioverbs::plotTotalReads
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotTotalReads(bcb_small)
NULL



#' @rdname plotTotalReads
#' @name plotTotalReads
#' @importFrom bioverbs plotTotalReads
#' @usage plotTotalReads(object, ...)
#' @export
NULL



plotTotalReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "total reads"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        p <- metrics(object) %>%
            ggplot(
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalReads") / 1e6L,
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
                y = "reads per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            # Convert limit to per million
            if (limit < 1e6L) {
                # nocov start
                warning("`limit`: Use absolute value, not per million")
                # nocov end
            } else {
                limit <- limit / 1e6L
            }
            if (limit > 1L) {
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



#' @rdname plotTotalReads
#' @export
setMethod(
    f = "plotTotalReads",
    signature = signature("bcbioRNASeq"),
    definition = plotTotalReads.bcbioRNASeq
)
