#' @name plotRRNAMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotRRNAMappingRate
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotRRNAMappingRate(bcb)
NULL



#' @importFrom basejump plotRRNAMappingRate
#' @aliases NULL
#' @export
basejump::plotRRNAMappingRate



plotRRNAMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.1,
        fill,
        flip,
        title = "rRNA mapping rate"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(object, interestingGroups)
        interestingGroups(object) <- interestingGroups
        assert_is_a_number(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsStringOrNULL(title)

        data <- metrics(object)

        # Warn and early return if rRNA rate was not calculated.
        if (!"rrnaRate" %in% colnames(data)) {
            warning("rRNA mapping rate was not calculated. Skipping plot.")
            return(invisible())
        }

        p <- ggplot(
            data = data,
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

formals(plotRRNAMappingRate.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotRRNAMappingRate.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    f = "plotRRNAMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotRRNAMappingRate.bcbioRNASeq
)
