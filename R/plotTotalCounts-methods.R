#' Plot Total Counts
#'
#' @name plotTotalCounts
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param normalized `boolean`. Use raw counts (`FALSE`) or DESeq2 normailzed
#'   counts (`TRUE`). Only applies to gene-level counts.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotTotalCounts(bcb_small)
NULL



#' @rdname plotTotalCounts
#' @export
setMethod(
    "plotTotalCounts",
    signature("bcbioRNASeq"),
    function(
        object,
        normalized = FALSE,
        interestingGroups = NULL,
        limit = 10e6L,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "total counts"
    ) {
        validObject(object)
        assert_is_a_bool(normalized)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        interestingGroups(object) <- interestingGroups
        assertIsAnImplicitInteger(limit)
        assert_all_are_non_negative(limit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        counts <- counts(object, normalized = normalized)
        data <- sampleData(object)
        data[["totalCounts"]] <- colSums(counts)

        p <- ggplot(
                data = as.data.frame(data),
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalCounts") / 1e6L,
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
                y = "counts per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            # Convert limit to per million.
            if (limit < 1e6L) {
                # nocov start
                warning("`limit`: Use absolute value, not per million")
                # nocov end
            } else {
                limit <- limit / 1e6L
            }
            if (limit > 1L) {
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
)
