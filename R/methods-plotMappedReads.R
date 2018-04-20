#' Plot Mapped Reads
#'
#' @name plotMappedReads
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotMappedReads(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotMappedReads
#' @export
setMethod(
    "plotMappedReads",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        passLimit = 20L,
        warnLimit = 10L,
        fill = scale_fill_hue(),
        flip = TRUE,
        title = "mapped reads"
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assertIsAnImplicitInteger(passLimit)
        assert_all_are_non_negative(passLimit)
        assertIsAnImplicitInteger(warnLimit)
        assert_all_are_non_negative(warnLimit)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        metrics <- metrics(object) %>%
            uniteInterestingGroups(interestingGroups)

        p <- ggplot(
            data = metrics,
            mapping = aes_(
                x = ~sampleName,
                y = ~mappedReads / 1e6L,
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
                y = "mapped reads per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (is_positive(passLimit)) {
            p <- p + .qcPassLine(passLimit)
        }
        if (is_positive(warnLimit)) {
            p <- p + .qcWarnLine(warnLimit)
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
