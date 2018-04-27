#' Plot Total Reads
#'
#' @name plotTotalReads
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotTotalReads(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotTotalReads
#' @export
setMethod(
    "plotTotalReads",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 10L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "total reads"
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
                y = ~totalReads / 1e6L,
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
                y = "reads per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

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
