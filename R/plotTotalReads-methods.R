#' @name plotTotalReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotTotalReads
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotTotalReads(bcb)
NULL



#' @importFrom basejump plotTotalReads
#' @aliases NULL
#' @export
basejump::plotTotalReads



plotTotalReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        perMillion = TRUE,
        fill = getOption("basejump.fill", NULL),
        flip = getOption("basejump", TRUE),
        title = "total reads"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assert_is_a_bool(perMillion)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsStringOrNULL(title)

        data <- metrics(object)

        # Convert to per million, if desired.
        yLab <- "reads"
        if (isTRUE(perMillion)) {
            data <- mutate(data, totalReads = !!sym("totalReads") / 1e6L)
            yLab <- paste(yLab, "per million")
        }

        p <- ggplot(
                data = data,
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalReads"),
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
                y = yLab,
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            if (isTRUE(perMillion)) {
                if (limit < 1e6L) {
                    warning("`limit`: Use absolute value, not per million.")
                } else {
                    limit <- limit / 1e6L
                }
            }
            p <- p + basejump_geom_abline(yintercept = limit)
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
