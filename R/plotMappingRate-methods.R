#' @name plotMappingRate
#' @inherit basejump::plotMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams params
#' @inheritParams basejump::params
#'
#' @examples
#' data(bcb)
#' plotMappingRate(bcb)
NULL



#' @importFrom basejump plotMappingRate
#' @aliases NULL
#' @export
basejump::plotMappingRate



# bcbioRNASeq ==================================================================
plotMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.7,
        fill = getOption("basejump.fill", NULL),
        flip = getOption("basejump", TRUE),
        title = "mapping rate"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assert_is_a_number(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / !!sym("totalReads") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            ylim(0L, 100L) +
            labs(
                title = title,
                x = NULL,
                y = "mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            # Convert to percentage
            if (limit > 1L) {
                # nocov start
                warning("`limit`: Use ratio (0-1) instead of percentage.")
                # nocov end
            } else {
                limit <- limit * 100L
            }
            if (limit < 100L) {
                p <- p + basejump_geom_abline(yintercept = limit)
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



#' @rdname plotMappingRate
#' @export
setMethod(
    f = "plotMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotMappingRate.bcbioRNASeq
)
