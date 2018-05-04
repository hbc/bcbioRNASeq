#' Plot 5'->3' Bias
#'
#' @name plot5Prime3PrimeBias
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plot5Prime3PrimeBias(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    "plot5Prime3PrimeBias",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        limit = 2L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "5'->3' bias"
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
            mapping = aes_string(
                x = "sampleName",
                y = "x5x3Bias",
                fill = "interestingGroups"
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = "5'->3' bias",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            p <- p + bcbio_geom_abline(yintercept = limit)
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
