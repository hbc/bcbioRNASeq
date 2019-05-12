#' @name plot5Prime3PrimeBias
#' @inherit bioverbs::plot5Prime3PrimeBias
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @details
#' It is common to see a small amount of bias, especially if polyA enrichment
#' was performed, or if there is any sample degradation. If a large amount of
#' bias is observed here, be sure to analyze the samples with a Bioanalyzer and
#' check the RIN scores.
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `ggplot`.
#'
#' @examples
#' plot5Prime3PrimeBias(bcb_small)
NULL



#' @rdname plot5Prime3PrimeBias
#' @name plot5Prime3PrimeBias
#' @importFrom bioverbs plot5Prime3PrimeBias
#' @usage plot5Prime3PrimeBias(object, ...)
#' @export
NULL



plot5Prime3PrimeBias.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups,
        limit = 0L,
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "5'->3' bias"
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
                y = !!sym("x5x3Bias"),
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
                y = "5'->3' bias",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(limit)) {
            p <- p + bcbio_geom_abline(yintercept = limit)  # nocov
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



#' @rdname plot5Prime3PrimeBias
#' @export
setMethod(
    f = "plot5Prime3PrimeBias",
    signature = signature("bcbioRNASeq"),
    definition = plot5Prime3PrimeBias.bcbioRNASeq
)
