#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotMappedReads
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotMappedReads(bcb)
NULL



#' @importFrom basejump plotMappedReads
#' @aliases NULL
#' @export
basejump::plotMappedReads



plotMappedReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        fill,
        flip,
        title = "mapped reads"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsStringOrNULL(title)

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / 1e6L,
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
                y = "mapped reads per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            # Convert limit to per million.
            if (limit < 1e6L) {
                # nocov start
                warning("`limit`: Use absolute value, not per million.")
                # nocov end
            } else {
                limit <- limit / 1e6L
            }
            if (limit >= 1L) {
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

formals(plotMappedReads.bcbioRNASeq)[["fill"]] <- fillDiscrete
formals(plotMappedReads.bcbioRNASeq)[["flip"]] <- flip



#' @rdname plotMappedReads
#' @export
setMethod(
    f = "plotMappedReads",
    signature = signature("bcbioRNASeq"),
    definition = plotMappedReads.bcbioRNASeq
)
