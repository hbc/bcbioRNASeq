#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @name plotRRNAMappingRate
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotRRNAMappingRate(bcb_small)
NULL



# Constructors =================================================================
.plotRRNAMappingRate <- function(
    object,
    interestingGroups,
    warnLimit = 10L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = "rRNA mapping rate"
) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
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
            y = ~rrnaRate * 100L,
            fill = ~interestingGroups
        )
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "rRNA mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        )

    if (is_positive(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
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



# Methods ======================================================================
#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeq"),
    .plotRRNAMappingRate
)
